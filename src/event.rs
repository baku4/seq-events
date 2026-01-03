/// Parsing event. References are valid until the next `next_event()` call.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event<'a> {
    /// Next record starting (signals end of previous record).
    /// Not emitted before the first record.
    NextRecord,
    /// Record ID chunk (may be partial if spanning buffer boundary).
    IdChunk(&'a [u8]),
    /// Sequence data chunk.
    SeqChunk(&'a [u8]),
    /// Quality scores chunk (FASTQ only).
    QualChunk(&'a [u8]),
}
