use std::io;
use thiserror::Error;

/// Errors from sequence parsing.
#[derive(Debug, Error)]
pub enum ReaderError {
    #[error("IO error: {0}")]
    Io(#[from] io::Error),

    #[error("Invalid format: {message}")]
    InvalidFormat { message: String },
}
