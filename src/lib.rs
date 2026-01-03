mod error;
mod event;
mod fasta;
mod fastq;

pub use error::EventSeqReaderError;
pub use event::Event;
pub use fasta::FastaReader;
pub use fastq::FastqReader;
