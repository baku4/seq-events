//! Zero-copy streaming FASTA/FASTQ parser.
//!
//! ```no_run
//! use seq_events::{FastaReader, Event};
//! use std::fs::File;
//!
//! let mut reader = FastaReader::new(File::open("seq.fasta").unwrap());
//! while let Some(Ok(event)) = reader.next_event() {
//!     match event {
//!         Event::StartRecord => {}
//!         Event::IdChunk(id) => {}
//!         Event::SeqChunk(seq) => {}
//!         Event::QualChunk(_) => unreachable!(),
//!     }
//! }
//! ```

mod error;
mod event;
mod fasta;
mod fastq;

pub use error::ReaderError;
pub use event::Event;
pub use fasta::FastaReader;
pub use fastq::FastqReader;
