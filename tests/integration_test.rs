use std::fs::File;
use std::path::Path;

use flate2::read::GzDecoder;
use seq_events::{Event, FastaReader, FastqReader};

const TEST_DATA_DIR: &str = "tests/test_data";

fn fasta_dir() -> std::path::PathBuf {
    Path::new(TEST_DATA_DIR).join("fasta")
}

fn fastq_dir() -> std::path::PathBuf {
    Path::new(TEST_DATA_DIR).join("fastq")
}

fn count_fasta_stats<R: std::io::Read>(mut reader: FastaReader<R>) -> (usize, usize, Vec<String>) {
    let mut record_count = 0;
    let mut total_seq_len = 0;
    let mut record_ids = Vec::new();
    let mut current_id = Vec::new();

    while let Some(event) = reader.next_event() {
        match event.expect("Failed to parse FASTA") {
            Event::StartRecord => {
                record_count += 1;
                current_id.clear();
            }
            Event::IdChunk(chunk) => {
                current_id.extend_from_slice(chunk);
            }
            Event::SeqChunk(bases) => {
                if record_ids.len() < record_count {
                    record_ids.push(String::from_utf8_lossy(&current_id).to_string());
                    current_id.clear();
                }
                total_seq_len += bases.len();
            }
            Event::QualChunk(_) => unreachable!(),
        }
    }

    if record_ids.len() < record_count {
        record_ids.push(String::from_utf8_lossy(&current_id).to_string());
    }

    (record_count, total_seq_len, record_ids)
}

fn count_fastq_stats<R: std::io::Read>(
    mut reader: FastqReader<R>,
) -> (usize, usize, usize, Vec<String>) {
    let mut record_count = 0;
    let mut total_seq_len = 0;
    let mut total_qual_len = 0;
    let mut record_ids = Vec::new();
    let mut current_id = Vec::new();

    while let Some(event) = reader.next_event() {
        match event.expect("Failed to parse FASTQ") {
            Event::StartRecord => {
                record_count += 1;
                current_id.clear();
            }
            Event::IdChunk(chunk) => {
                current_id.extend_from_slice(chunk);
            }
            Event::SeqChunk(bases) => {
                if record_ids.len() < record_count {
                    record_ids.push(String::from_utf8_lossy(&current_id).to_string());
                    current_id.clear();
                }
                total_seq_len += bases.len();
            }
            Event::QualChunk(quals) => {
                total_qual_len += quals.len();
            }
        }
    }

    if record_ids.len() < record_count {
        record_ids.push(String::from_utf8_lossy(&current_id).to_string());
    }

    (record_count, total_seq_len, total_qual_len, record_ids)
}

#[test]
fn test_fasta_influenza_lf() {
    let path = fasta_dir().join("influenza.fasta");
    if !path.exists() {
        panic!("File not found: {}", path.display());
    }

    let file = File::open(&path).unwrap();
    let (record_count, total_seq_len, record_ids) = count_fasta_stats(FastaReader::new(file));

    assert_eq!(record_count, 8);
    assert_eq!(total_seq_len, 13627);
    assert!(record_ids[0].starts_with("NC_007373.1"));
    assert!(record_ids[7].starts_with("NC_007370.1"));
}

#[test]
fn test_fasta_influenza_crlf() {
    let path = fasta_dir().join("influenza_crlf.fasta");
    if !path.exists() {
        panic!("File not found: {}", path.display());
    }

    let file = File::open(&path).unwrap();
    let (record_count, total_seq_len, record_ids) = count_fasta_stats(FastaReader::new(file));

    assert_eq!(record_count, 8);
    assert_eq!(total_seq_len, 13627);
    assert!(record_ids[0].starts_with("NC_007373.1"));
    assert!(record_ids[7].starts_with("NC_007370.1"));
}

#[test]
fn test_fasta_human_chr22_gzip() {
    let path = fasta_dir().join("human_chr22.fasta.gz");
    if !path.exists() {
        panic!("File not found: {}", path.display());
    }

    let file = File::open(&path).unwrap();
    let decoder = GzDecoder::new(file);
    let (record_count, total_seq_len, record_ids) = count_fasta_stats(FastaReader::new(decoder));

    assert_eq!(record_count, 1);
    assert!(total_seq_len > 40_000_000);
    assert!(record_ids[0].contains("22") || record_ids[0].contains("NC_"));
}

#[test]
fn test_fastq_sample_lf() {
    let path = fastq_dir().join("sample.fastq");
    if !path.exists() {
        panic!("File not found: {}", path.display());
    }

    let file = File::open(&path).unwrap();
    let (record_count, total_seq_len, total_qual_len, record_ids) =
        count_fastq_stats(FastqReader::new(file));

    assert_eq!(record_count, 5);
    assert_eq!(total_seq_len, 252);
    assert_eq!(total_seq_len, total_qual_len);
    assert!(record_ids[0].starts_with("read1"));
    assert!(record_ids[4].starts_with("read5"));
}

#[test]
fn test_fastq_sample_crlf() {
    let path = fastq_dir().join("sample_crlf.fastq");
    if !path.exists() {
        panic!("File not found: {}", path.display());
    }

    let file = File::open(&path).unwrap();
    let (record_count, total_seq_len, total_qual_len, record_ids) =
        count_fastq_stats(FastqReader::new(file));

    assert_eq!(record_count, 5);
    assert_eq!(total_seq_len, 252);
    assert_eq!(total_seq_len, total_qual_len);
    assert!(record_ids[0].starts_with("read1"));
    assert!(record_ids[4].starts_with("read5"));
}

#[test]
fn test_fastq_sample_gzip() {
    let path = fastq_dir().join("sample.fastq.gz");
    if !path.exists() {
        panic!("File not found: {}", path.display());
    }

    let file = File::open(&path).unwrap();
    let decoder = GzDecoder::new(file);
    let (record_count, total_seq_len, total_qual_len, record_ids) =
        count_fastq_stats(FastqReader::new(decoder));

    assert_eq!(record_count, 5);
    assert_eq!(total_seq_len, 252);
    assert_eq!(total_seq_len, total_qual_len);
    assert!(record_ids[0].starts_with("read1"));
    assert!(record_ids[4].starts_with("read5"));
}

#[test]
fn test_fasta_lf_vs_crlf_consistency() {
    let lf_path = fasta_dir().join("influenza.fasta");
    let crlf_path = fasta_dir().join("influenza_crlf.fasta");

    if !lf_path.exists() {
        panic!("File not found: {}", lf_path.display());
    }
    if !crlf_path.exists() {
        panic!("File not found: {}", crlf_path.display());
    }

    let lf_stats = count_fasta_stats(FastaReader::new(File::open(&lf_path).unwrap()));
    let crlf_stats = count_fasta_stats(FastaReader::new(File::open(&crlf_path).unwrap()));

    assert_eq!(lf_stats, crlf_stats);
}

#[test]
fn test_fastq_lf_vs_crlf_consistency() {
    let lf_path = fastq_dir().join("sample.fastq");
    let crlf_path = fastq_dir().join("sample_crlf.fastq");

    if !lf_path.exists() {
        panic!("File not found: {}", lf_path.display());
    }
    if !crlf_path.exists() {
        panic!("File not found: {}", crlf_path.display());
    }

    let lf_stats = count_fastq_stats(FastqReader::new(File::open(&lf_path).unwrap()));
    let crlf_stats = count_fastq_stats(FastqReader::new(File::open(&crlf_path).unwrap()));

    assert_eq!(lf_stats, crlf_stats);
}
