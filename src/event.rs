/// Parsing event. References are valid until the next `next_event()` call.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event<'a> {
    /// New record started (implicitly ends previous record).
    StartRecord,
    /// Record ID chunk (may be partial if spanning buffer boundary).
    IdChunk(&'a [u8]),
    /// Sequence data chunk.
    SeqChunk(&'a [u8]),
    /// Quality scores chunk (FASTQ only).
    QualChunk(&'a [u8]),
}
