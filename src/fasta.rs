use std::io::{BufRead, BufReader, Read};

use memchr::{memchr, memchr3};

use crate::error::EventSeqReaderError;
use crate::event::Event;

const DEFAULT_BUFFER_SIZE: usize = 128 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Start,
    Id,
    Sequence,
}

pub struct FastaReader<R> {
    reader: BufReader<R>,
    pending_consume: usize,
    state: State,
}

impl<R: Read> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        Self::with_capacity(DEFAULT_BUFFER_SIZE, reader)
    }

    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            pending_consume: 0,
            state: State::Start,
        }
    }

    pub fn next_event(&mut self) -> Option<Result<Event<'_>, EventSeqReaderError>> {
        loop {
            if self.pending_consume > 0 {
                self.reader.consume(self.pending_consume);
                self.pending_consume = 0;
            }

            let buf = match self.reader.fill_buf() {
                Ok(b) if b.is_empty() => return None,
                Ok(b) => b,
                Err(e) => return Some(Err(e.into())),
            };

            let buf_ptr = buf.as_ptr();
            let buf_len = buf.len();

            match self.state {
                State::Start => {
                    let first_non_ws = buf.iter().position(|&b| b != b'\n' && b != b'\r');

                    match first_non_ws {
                        Some(0) => {
                            if buf[0] == b'>' {
                                self.state = State::Id;
                                self.pending_consume = 1;
                                return Some(Ok(Event::StartRecord));
                            } else {
                                return Some(Err(EventSeqReaderError::InvalidFormat {
                                    message: format!(
                                        "Expected '>' at start of FASTA record, found '{}'",
                                        buf[0] as char
                                    ),
                                }));
                            }
                        }
                        Some(pos) => {
                            self.pending_consume = pos;
                            continue;
                        }
                        None => {
                            self.pending_consume = buf_len;
                            continue;
                        }
                    }
                }

                State::Id => {
                    if let Some(newline_pos) = memchr(b'\n', buf) {
                        let end = if newline_pos > 0 && buf[newline_pos - 1] == b'\r' {
                            newline_pos - 1
                        } else {
                            newline_pos
                        };

                        self.state = State::Sequence;
                        self.pending_consume = newline_pos + 1;

                        if end > 0 {
                            let slice = unsafe { std::slice::from_raw_parts(buf_ptr, end) };
                            return Some(Ok(Event::IdChunk(slice)));
                        } else {
                            continue;
                        }
                    } else {
                        self.pending_consume = buf_len;
                        let slice = unsafe { std::slice::from_raw_parts(buf_ptr, buf_len) };
                        return Some(Ok(Event::IdChunk(slice)));
                    }
                }

                State::Sequence => {
                    let first_byte = buf[0];

                    if first_byte == b'\n' {
                        self.pending_consume = 1;
                        continue;
                    }
                    if first_byte == b'\r' {
                        self.pending_consume = if buf_len > 1 && buf[1] == b'\n' { 2 } else { 1 };
                        continue;
                    }
                    if first_byte == b'>' {
                        self.state = State::Id;
                        self.pending_consume = 1;
                        return Some(Ok(Event::StartRecord));
                    }

                    let chunk_end = memchr3(b'\n', b'\r', b'>', buf).unwrap_or(buf_len);
                    if chunk_end == 0 {
                        self.pending_consume = 1;
                        continue;
                    }

                    self.pending_consume = chunk_end;
                    let slice = unsafe { std::slice::from_raw_parts(buf_ptr, chunk_end) };
                    return Some(Ok(Event::SeqChunk(slice)));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_single_record() {
        let data = b">seq1 description\nACGT\nTGCA\n";
        let mut reader = FastaReader::new(Cursor::new(&data[..]));

        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::StartRecord));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::IdChunk(id) if id == b"seq1 description"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"ACGT"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"TGCA"));
        assert!(reader.next_event().is_none());
    }

    #[test]
    fn test_multiple_records() {
        let data = b">seq1\nACGT\n>seq2\nTGCA\n";
        let mut reader = FastaReader::new(Cursor::new(&data[..]));

        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::StartRecord));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::IdChunk(id) if id == b"seq1"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"ACGT"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::StartRecord));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::IdChunk(id) if id == b"seq2"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"TGCA"));
        assert!(reader.next_event().is_none());
    }

    #[test]
    fn test_crlf_line_endings() {
        let data = b">seq1\r\nACGT\r\nTGCA\r\n";
        let mut reader = FastaReader::new(Cursor::new(&data[..]));

        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::StartRecord));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::IdChunk(id) if id == b"seq1"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"ACGT"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"TGCA"));
        assert!(reader.next_event().is_none());
    }

    #[test]
    fn test_small_buffer() {
        let data = b">seq1\nACGTACGTACGT\n";
        let mut reader = FastaReader::with_capacity(4, Cursor::new(&data[..]));

        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::StartRecord));

        let mut id = Vec::new();
        let mut seq = Vec::new();

        loop {
            match reader.next_event() {
                Some(Ok(Event::IdChunk(chunk))) => id.extend_from_slice(chunk),
                Some(Ok(Event::SeqChunk(chunk))) => seq.extend_from_slice(chunk),
                Some(Ok(_)) => panic!("Unexpected event"),
                Some(Err(e)) => panic!("Error: {}", e),
                None => break,
            }
        }

        assert_eq!(&id, b"seq1");
        assert_eq!(&seq, b"ACGTACGTACGT");
    }
}
