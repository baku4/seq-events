#![doc = include_str!("../README.md")]

mod error;
mod event;
mod fasta;
mod fastq;

pub use error::ReaderError;
pub use event::Event;
pub use fasta::FastaReader;
pub use fastq::FastqReader;
