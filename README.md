# seq-events

A minimal, zero-copy streaming parser for FASTA/FASTQ files.

Emits SAX-style events as it reads, never buffering entire records. Works with any `Read` source including gzip streams.

## Usage

```rust,no_run
use seq_events::{FastaReader, Event};
use std::fs::File;

let file = File::open("sequences.fasta").unwrap();
let mut reader = FastaReader::new(file);

while let Some(Ok(event)) = reader.next_event() {
    match event {
        Event::NextRecord => {}
        Event::IdChunk(id) => println!("{}", String::from_utf8_lossy(id)),
        Event::SeqChunk(seq) => println!("{} bp", seq.len()),
        Event::QualChunk(_) => {}
    }
}
```

```rust,no_run
use seq_events::{FastqReader, Event};
use std::fs::File;

let file = File::open("reads.fastq").unwrap();
let mut reader = FastqReader::new(file);

while let Some(Ok(event)) = reader.next_event() {
    match event {
        Event::NextRecord => {}
        Event::IdChunk(id) => println!("{}", String::from_utf8_lossy(id)),
        Event::SeqChunk(seq) => println!("{} bp", seq.len()),
        Event::QualChunk(qual) => println!("{} qual", qual.len()),
    }
}
```

## Events

- `NextRecord` - Emitted between records (not before the first)
- `IdChunk(&[u8])` - Record identifier (may span multiple chunks)
- `SeqChunk(&[u8])` - Sequence bases
- `QualChunk(&[u8])` - Phred quality scores, ASCII-encoded (FASTQ only)
