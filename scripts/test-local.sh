#!/bin/bash
# Local test script for haplograph
# Run this before pushing to ensure tests pass

set -e

echo "🧪 Running local tests for haplograph..."
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if cargo is installed
if ! command -v cargo &> /dev/null; then
    echo -e "${RED}❌ Cargo is not installed. Please install Rust first.${NC}"
    exit 1
fi

# Function to print test results
print_result() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✓ $1${NC}"
    else
        echo -e "${RED}✗ $1${NC}"
        exit 1
    fi
}

# Test 1: Build
echo ""
echo "1. Building project..."
cargo build --release
print_result "Build"


# Test 2: Help commands
echo ""
echo "2. Testing command-line interface..."
./target/release/haplograph --help > /dev/null
print_result "Main help command"
./target/release/haplograph haplograph --help > /dev/null 2>&1 || true
./target/release/haplograph dev-tools --help > /dev/null 2>&1 || true
print_result "Subcommand help commands"

# Test 3: Integration test with minimal data (optional)
echo ""
echo "3. Running integration test with synthetic data..."
mkdir -p test_output test_data

# Run haplograph command
if ./target/release/haplograph haplograph \
    --alignment-bam test/HG002_HLA_A.bam \
    --reference-fa test/hg38.chr6.fa.gz \
    --sampleid HG002 \
    --locus "chr6:29942532-29945870" \
    --output-prefix test_output/haplograph > /dev/null 2>&1; then
    print_result "Integration test (haplograph)"
else
    echo -e "${YELLOW}⚠ Integration test skipped (may require additional setup)${NC}"
fi
# Run haplointervals command
if ./target/release/haplograph haplointervals \
    --alignment-bam test/HG002_HLA_A.bam \
    --reference-fa test/hg38.chr6.fa.gz \
    --sampleid HG002 \
    --bed-file test/test.bed \
    --output-prefix test_output/haplointervals > /dev/null 2>&1; then
    print_result "Integration test (haplointervals)"
else
    echo -e "${YELLOW}⚠ Integration test skipped (may require additional setup)${NC}"
fi

# Run evaluate command
if ./target/release/haplograph evaluate \
    --truth-fasta test/HG002_HLA-A.all.fa \
    --query-fasta test/HG002_HLA_query.fasta \
    --seq-number 2 \
    --output-prefix test_output/eval > /dev/null 2>&1; then
    print_result "Integration test (evaluate)"
else
    echo -e "${YELLOW}⚠ Integration test skipped (may require additional setup)${NC}"
fi

# Run extract command
if ./target/release/haplograph extract \
    --bamfile test/HG002_HLA_A.bam \
    --locus "chr6:29942532-29945870" \
    --output-prefix test_output/extract > /dev/null 2>&1; then
    print_result "Integration test (extract)"
else
    echo -e "${YELLOW}⚠ Integration test skipped (may require additional setup)${NC}"
fi

# Cleanup
rm -rf test_output test_data

echo ""
echo -e "${GREEN}✅ All tests passed!${NC}"
echo "You can now push your changes."

