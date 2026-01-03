use std::io::{BufRead, BufReader, Read};

use memchr::{memchr, memchr2, memchr3};

use crate::error::ReaderError;
use crate::event::Event;

const DEFAULT_BUFFER_SIZE: usize = 128 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    Start,
    Id,
    Sequence,
    Plus,
    Quality,
}

/// Zero-copy streaming FASTQ parser.
pub struct FastqReader<R> {
    reader: BufReader<R>,
    pending_consume: usize,
    state: State,
    seq_len: usize,
    qual_len: usize,
    first_record: bool,
}

impl<R: Read> FastqReader<R> {
    /// Creates a reader with default 128 KiB buffer.
    pub fn new(reader: R) -> Self {
        Self::with_capacity(DEFAULT_BUFFER_SIZE, reader)
    }

    /// Creates a reader with specified buffer capacity.
    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Self {
            reader: BufReader::with_capacity(capacity, reader),
            pending_consume: 0,
            state: State::Start,
            seq_len: 0,
            qual_len: 0,
            first_record: true,
        }
    }

    /// Returns the next event, or `None` at EOF.
    pub fn next_event(&mut self) -> Option<Result<Event<'_>, ReaderError>> {
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
                            if buf[0] == b'@' {
                                let is_first = self.first_record;
                                self.first_record = false;
                                self.state = State::Id;
                                self.pending_consume = 1;
                                self.seq_len = 0;
                                self.qual_len = 0;
                                if is_first {
                                    continue; // First record - no event
                                }
                                return Some(Ok(Event::NextRecord));
                            } else {
                                return Some(Err(ReaderError::InvalidFormat {
                                    message: format!(
                                        "Expected '@' at start of FASTQ record, found '{}'",
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
                    if buf[0] == b'+' {
                        if let Some(newline_pos) = memchr(b'\n', buf) {
                            self.pending_consume = newline_pos + 1;
                        } else {
                            self.pending_consume = buf_len;
                        }
                        self.state = State::Plus;
                        continue;
                    }

                    if buf[0] == b'\n' {
                        self.pending_consume = 1;
                        continue;
                    }
                    if buf[0] == b'\r' {
                        self.pending_consume = if buf_len > 1 && buf[1] == b'\n' { 2 } else { 1 };
                        continue;
                    }

                    let chunk_end = memchr3(b'\n', b'\r', b'+', buf).unwrap_or(buf_len);

                    if chunk_end == 0 {
                        self.pending_consume = 1;
                        continue;
                    }

                    self.pending_consume = chunk_end;
                    self.seq_len += chunk_end;
                    let slice = unsafe { std::slice::from_raw_parts(buf_ptr, chunk_end) };
                    return Some(Ok(Event::SeqChunk(slice)));
                }

                State::Plus => {
                    if buf[0] == b'\n' {
                        self.pending_consume = 1;
                        self.state = State::Quality;
                        continue;
                    }
                    if buf[0] == b'\r' {
                        self.pending_consume = if buf_len > 1 && buf[1] == b'\n' { 2 } else { 1 };
                        self.state = State::Quality;
                        continue;
                    }
                    self.state = State::Quality;
                    continue;
                }

                State::Quality => {
                    if buf[0] == b'\n' {
                        self.pending_consume = 1;
                        continue;
                    }
                    if buf[0] == b'\r' {
                        self.pending_consume = if buf_len > 1 && buf[1] == b'\n' { 2 } else { 1 };
                        continue;
                    }

                    let remaining = self.seq_len.saturating_sub(self.qual_len);

                    if remaining == 0 {
                        self.state = State::Start;
                        continue;
                    }

                    let newline_pos = memchr2(b'\n', b'\r', buf);
                    let chunk_end = match newline_pos {
                        Some(pos) => pos.min(remaining),
                        None => buf_len.min(remaining),
                    };

                    if chunk_end == 0 {
                        self.pending_consume = 1;
                        continue;
                    }

                    self.pending_consume = chunk_end;
                    self.qual_len += chunk_end;
                    let slice = unsafe { std::slice::from_raw_parts(buf_ptr, chunk_end) };

                    if self.qual_len >= self.seq_len {
                        self.state = State::Start;
                    }

                    return Some(Ok(Event::QualChunk(slice)));
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
        let data = b"@read1 description\nACGT\n+\nIIII\n";
        let mut reader = FastqReader::new(Cursor::new(&data[..]));

        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::IdChunk(id) if id == b"read1 description"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"ACGT"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::QualChunk(q) if q == b"IIII"));
        assert!(reader.next_event().is_none());
    }

    #[test]
    fn test_multiple_records() {
        let data = b"@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nHHHH\n";
        let mut reader = FastqReader::new(Cursor::new(&data[..]));

        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::IdChunk(id) if id == b"read1"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"ACGT"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::QualChunk(q) if q == b"IIII"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::NextRecord));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::IdChunk(id) if id == b"read2"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"TGCA"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::QualChunk(q) if q == b"HHHH"));
        assert!(reader.next_event().is_none());
    }

    #[test]
    fn test_crlf_line_endings() {
        let data = b"@read1\r\nACGT\r\n+\r\nIIII\r\n";
        let mut reader = FastqReader::new(Cursor::new(&data[..]));

        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::IdChunk(id) if id == b"read1"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::SeqChunk(s) if s == b"ACGT"));
        assert!(matches!(reader.next_event().unwrap().unwrap(), Event::QualChunk(q) if q == b"IIII"));
        assert!(reader.next_event().is_none());
    }

    #[test]
    fn test_small_buffer() {
        let data = b"@read1\nACGTACGT\n+\nIIIIIIII\n";
        let mut reader = FastqReader::with_capacity(4, Cursor::new(&data[..]));

        let mut id = Vec::new();
        let mut seq = Vec::new();
        let mut qual = Vec::new();

        loop {
            match reader.next_event() {
                Some(Ok(Event::IdChunk(chunk))) => id.extend_from_slice(chunk),
                Some(Ok(Event::SeqChunk(chunk))) => seq.extend_from_slice(chunk),
                Some(Ok(Event::QualChunk(chunk))) => qual.extend_from_slice(chunk),
                Some(Ok(_)) => panic!("Unexpected event"),
                Some(Err(e)) => panic!("Error: {}", e),
                None => break,
            }
        }

        assert_eq!(&id, b"read1");
        assert_eq!(&seq, b"ACGTACGT");
        assert_eq!(&qual, b"IIIIIIII");
    }
}
