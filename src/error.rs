use std::io;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum EventSeqReaderError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),

    #[error("Invalid format: {message}")]
    InvalidFormat { message: String },
}
