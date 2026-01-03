#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEST_DATA_DIR="$SCRIPT_DIR/test_data"
FASTA_DIR="$TEST_DATA_DIR/fasta"
FASTQ_DIR="$TEST_DATA_DIR/fastq"

mkdir -p "$FASTA_DIR" "$FASTQ_DIR"

echo "=== Downloading test data ==="

echo "Downloading Influenza A genome..."
INFLUENZA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/865/085/GCF_000865085.1_ViralMultiSegProj15622/GCF_000865085.1_ViralMultiSegProj15622_genomic.fna.gz"
curl -sL "$INFLUENZA_URL" | gunzip > "$FASTA_DIR/influenza.fasta"
echo "  -> influenza.fasta ($(wc -l < "$FASTA_DIR/influenza.fasta") lines)"

echo "Creating CRLF version..."
sed 's/$/\r/' "$FASTA_DIR/influenza.fasta" > "$FASTA_DIR/influenza_crlf.fasta"
echo "  -> influenza_crlf.fasta"

echo "Downloading Human chromosome 22..."
CHR22_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr22.fna.gz"
curl -sL "$CHR22_URL" -o "$FASTA_DIR/human_chr22.fasta.gz"
echo "  -> human_chr22.fasta.gz"

echo "Creating synthetic FASTQ sample..."
cat > "$FASTQ_DIR/sample.fastq" << 'EOF'
@read1 length=50
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2 length=50
TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATG
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@read3 length=50
AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAA
+
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
@read4 length=50
ATATATATATATATATATATATATATATATATATATATATATATATATATAT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@read5 length=50
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
+
EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
EOF
echo "  -> sample.fastq (5 reads)"

sed 's/$/\r/' "$FASTQ_DIR/sample.fastq" > "$FASTQ_DIR/sample_crlf.fastq"
echo "  -> sample_crlf.fastq"

gzip -k "$FASTQ_DIR/sample.fastq"
echo "  -> sample.fastq.gz"

echo ""
echo "=== Download complete ==="
echo "Test data location: $TEST_DATA_DIR"
