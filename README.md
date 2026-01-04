# seq-events

Zero-copy, event-driven streaming FASTA/FASTQ parser.

## Example

```rust
use seq_events::{FastaReader, Event};
use std::fs::File;

let file = File::open("sequences.fasta")?;
let mut reader = FastaReader::new(file);

while let Some(Ok(event)) = reader.next_event() {
    match event {
        Event::NextRecord => println!("---"),
        Event::IdChunk(id) => println!("ID: {}", String::from_utf8_lossy(id)),
        Event::SeqChunk(seq) => println!("Seq: {} bp", seq.len()),
        Event::QualChunk(_) => unreachable!(),
    }
}
```

For FASTQ:

```rust
use seq_events::{FastqReader, Event};

let mut reader = FastqReader::new(file);

while let Some(Ok(event)) = reader.next_event() {
    match event {
        Event::NextRecord => { /* new record */ }
        Event::IdChunk(id) => { /* record ID */ }
        Event::SeqChunk(seq) => { /* sequence bases */ }
        Event::QualChunk(qual) => { /* Phred scores, ASCII-encoded */ }
    }
}
```
