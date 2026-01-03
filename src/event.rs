#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Event<'a> {
    StartRecord,
    IdChunk(&'a [u8]),
    SeqChunk(&'a [u8]),
    QualChunk(&'a [u8]),
}
