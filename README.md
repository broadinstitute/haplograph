# Haplograph

A  bioinformatics toolkit for pangenome-guided local diploid assembly, methylation analysis and variant calling from lrWGS BAM files using graph-based approaches.

## Overview

Haplograph is designed for accurate haplotype and meplotype construction and variant calling in local genomic regions, particularly useful for complex region such as HLA/MHC analysis, and other highly polymorphic regions. The tool uses graph-based assembly to reconstruct haplotypes (methylation aware) and call variants with high accuracy.

## Features

- **Graph-based Haplotype Assembly**: Constructs sequence graphs from aligned reads
- **Variant Calling**: Generates VCF files with phased variants and methylation information
- **Methylation Analysis**: Incorporates methylation signals for improved haplotype resolution
- **Multiple Output Formats**: Supports GFA (Graph Fragment Assembly), VCF, and FASTA formats
- **Germline Haplotype Filtering**: Focus on major haplotypes for cleaner results
- **Read-based Phasing**: Uses read overlap information for accurate phasing
- **Comprehensive Evaluation**: Compare results against truth sets with QV scoring
- **Flexible Parameters**: Configurable window sizes, read support thresholds, and filtering options
- **Multiple Analysis Modes**: Support for continuous regions (haplograph) and interval-based analysis (haplointervals)
- **CI/CD Integration**: Automated testing and building with GitHub Actions

## Installation

### Prerequisites

- Rust 1.70+ (install from [rustup.rs](https://rustup.rs/))
- HTSlib development libraries
- Standard bioinformatics tools (samtools, minimap2) for file processing

### Build from Source

```bash
# Clone the repository
git clone https://github.com/broadinstitute/haplograph.git
cd haplograph

# Build the project
cargo build --release

# Install globally (optional)
cargo install --path .
```

## Usage

Haplograph provides several commands for different stages of haplotype analysis:

### Quick-start
```bash
haplograph haplograph \
    --alignment-bam input.bam \
    --reference-fa reference.fa \
    --sampleid SAMPLE001 \
    --locus chr6:29943661-29943700 \
    --output-prefix output/HLA_A
```
### 3. HaploPan (Pangenome-guided Analysis)

HaploPan leverages a pangenome FASTA to resolve complex haplotypes (e.g., duplications, deletions, and translocations) and then realigns reads to the best matching pangenome paths before running local assembly.

```bash
haplograph haplopan \
    --pangenome-fasta pangenome.fa \
    --alignment-bam input.bam \
    --sample-id SAMPLE001 \
    --rollingkmer-list 31,51,71 \
    --data-technology hifi \
    --window-size 100 \
    --output-prefix output/HLA_DRB1
```

**Parameters:**
- `--pangenome-fasta`: Pangenome FASTA with allele sequences
- `--alignment-bam`: Input BAM file
- `--sample-id`: Sample identifier
- `--rollingkmer-list`: Comma-separated k-mer sizes (default: 31)
- `--data-technology`: Sequencing technology: hifi, nanopore, or sr (default: hifi)
- `--window-size`: Window size for downstream graph/assembly (default: 100)
- `--output-prefix`: Output prefix for intermediate and final files

**Outputs:**
- `{prefix}.tmp.ref.fasta`: Best-matching pangenome allele sequences
- `{prefix}.tmp.bam` / `{prefix}.tmp.sorted.bam`: Realigned reads against the best alleles
- `{prefix}.final.fasta`: Final assembled haplotype sequences


### 2. Haplograph (Bam Continuous Region Analysis)

For analyzing continuous genomic regions (typically > 1kb):

Extract haplotypes from BAM files and build sequence graphs:

```bash
# Basic haplotype extraction
haplograph haplograph \
    --alignment-bam input.bam \
    --reference-fa reference.fa \
    --sampleid SAMPLE001 \
    --locus chr6:29943661-29943700 \
    --output-prefix output/HLA_A

# With custom parameters
haplograph haplograph \
    --alignment-bam input.bam \
    --reference-fa reference.fa \
    --sampleid SAMPLE001 \
    --locus chr6:29943661-29943700 \
    --frequency-min 0.05 \
    --min-reads 3 \
    --window-size 500 \
    --primary-only \
    --default-file-format gfa \
    --output-prefix output/HLA_A \
    --verbose
```

**Parameters:**
- `--alignment-bam`: Input BAM file
- `--reference-fa`: Reference FASTA file
- `--sampleid`: Sample identifier
- `--locus`: Genomic region (format: chr:start-end)
- `--var-frequency-min`: Minimum variant allele frequency (default: 0.01)
- `--min-reads`: Minimum supporting reads (default: 2)
- `--window-size`: Window size for analysis (default: 100)
- `--primary-only`: Use only primary alignments
- `--file-format`: Output format (fasta, gfa, or vcf, default: gfa)
- `--number-of-haplotypes`: Number of haplotypes to extract (default: 2)
- `--coverage-fold-threshold`: Heterozygous coverage fold threshold (default: 3.0)
- `--threshold-methyl-likelihood`: Methylation likelihood threshold (default: 0.5)
- `--detection-technology`: Sequencing technology: hifi or nanopore (default: hifi)

### 1. Haplointervals (Interval-based Analysis)

For analyzing multiple genomic intervals from a BED file (typically < 1kb each):

```bash
haplograph haplointervals \
    --alignment-bam input.bam \
    --reference-fa reference.fa \
    --sampleid SAMPLE001 \
    --bed-file intervals.bed \
    --output-prefix output/haplointervals
```

**Parameters:**
- `--alignment-bam`: Input BAM file
- `--reference-fa`: Reference FASTA file
- `--sampleid`: Sample identifier
- `--bed-file`: BED file containing genomic regions
- `--var-frequency-min`: Minimum variant allele frequency (default: 0.01)
- `--min-reads`: Minimum supporting reads (default: 2)
- `--window-size`: Maximum window size (default: 1000)
- `--primary-only`: Use only primary alignments
- `--threshold-methyl-likelihood`: Methylation likelihood threshold (default: 0.5)
- `--detection-technology`: Sequencing technology: hifi or nanopore (default: hifi)


### Helper Functions 1. Extract Sequences

Extract all sequences from a BAM file for a specific region:

```bash
haplograph extract \
    --bamfile input.bam \
    --locus chr6:29943661-29943700 \
    --output-prefix output/extracted
```

**Parameters:**
- `--bamfile`: Input BAM file
- `--locus`: Genomic region (format: chr:start-end)
- `--output-prefix`: Output prefix for FASTA file

### 6. Helper Functions 1. Evaluate Sequences, repurposed for genotyping

Evaluate haplotype accuracy against truth sets using QV (Quality Value) scoring:

```bash
haplograph evaluate \
    --truth-fasta truth_haplotypes.fasta \
    --query-fasta output/HLA_A_asm.fasta \
    --seq-number 2 \
    --output-prefix output/evaluation \
    --verbose
```

**Parameters:**
- `--truth-fasta`: Truth haplotype sequences (FASTA file)
- `--query-fasta`: Query haplotype sequences (FASTA file)
- `--seq-number`: Number of sequences to compare (default: 2)
- `--output-prefix`: Output prefix for evaluation TSV file


### 4. Haplotype Assembly (Dev Tools, implemented in Haplopan/ Haplograph / HaploInterval)

Assemble haplotypes from GFA files:

```bash
# Basic assembly
haplograph dev-tools assemble \
    --graph-gfa output/HLA_A.gfa \
    --output-prefix output/HLA_A_asm

# Germline-only assembly with specific haplotype count
haplograph dev-tools assemble \
    --graph-gfa output/HLA_A.gfa \
    --major-haplotype-only \
    --number-of-haplotypes 2 \
    --output-prefix output/HLA_A_asm \
    --verbose
```

**Parameters:**
- `--graph-gfa`: Input GFA file
- `--output-prefix`: Output prefix for assembled FASTA
- `--major-haplotype-only`: Focus on major haplotypes only
- `--number-of-haplotypes`: Number of haplotypes to extract (default: 2)
- `--fold-threshold`: Heterozygous coverage fold threshold (default: 3.0)

### 5. Variant Calling (Dev Tools, implemented in Haplopan/ Haplograph / HaploInterval)

Call variants from assembled haplotypes:

```bash
# Basic variant calling
haplograph dev-tools call \
    --gfa-file output/HLA_A_asm.gfa \
    --sampleid SAMPLE001 \
    --reference-fa reference.fa \
    --output-prefix output/HLA_A_variants

# With specific parameters
haplograph dev-tools call \
    --gfa-file output/HLA_A_asm.gfa \
    --sampleid SAMPLE001 \
    --reference-fa reference.fa \
    --maximum-haplotypes 2 \
    --fold-threshold 3.0 \
    --detection-technology hifi \
    --output-prefix output/HLA_A_variants \
    --verbose
```

**Parameters:**
- `--gfa-file`: Input GFA file from assembly
- `--sampleid`: Sample identifier
- `--reference-fa`: Reference FASTA file
- `--output-prefix`: Output prefix for VCF file
- `--maximum-haplotypes`: Maximum number of haplotypes (default: 2)
- `--fold-threshold`: Heterozygous coverage fold threshold (default: 3.0)
- `--detection-technology`: Sequencing technology: hifi or nanopore (default: hifi)


## Complete Workflow Example

```bash
# 1. Extract haplotypes and build graph
haplograph haplograph \
    --alignment-bam sample.bam \
    --reference-fa hg38.fa \
    --sampleid HG002 \
    --locus chr6:29943661-29943700 \
    --output-prefix output/HLA_A

# 2. Assemble haplotypes (optional, if using dev-tools)
haplograph dev-tools assemble \
    --graph-gfa output/HLA_A.gfa \
    --major-haplotype-only \
    --number-of-haplotypes 2 \
    --output-prefix output/HLA_A_asm

# 3. Call variants (optional, if using dev-tools)
haplograph dev-tools call \
    --gfa-file output/HLA_A_asm.gfa \
    --sampleid HG002 \
    --reference-fa hg38.fa \
    --output-prefix output/HLA_A_variants

# 4. Evaluate results (if truth available)
haplograph evaluate \
    --truth-fasta truth_HLA_A.fasta \
    --query-fasta output/HLA_A_asm.fasta \
    --seq-number 2 \
    --output-prefix output/evaluation
```

## Output Files

### Haplograph/Haplointervals
- `{prefix}.gfa`: Sequence graph in GFA format (default)
- `{prefix}.fasta`: Haplotype sequences (if fasta format selected)
- `{prefix}.vcf.gz`: Variant calls (if vcf format selected)
- `{prefix}.vcf.gz.tbi`: Tabix index for VCF file

### Assembly (dev-tools)
- `{prefix}.fasta`: Assembled haplotype sequences
- Contains multiple haplotype sequences with support information

### Variant Calling (dev-tools)
- `{prefix}.vcf.gz`: Compressed VCF file with called variants
- `{prefix}.vcf.gz.tbi`: Tabix index for the VCF file
- Includes phased variants with quality scores, allele depths (AD), variant allele frequencies (VAF), and methylation scores (MOD)

### Extract
- `{prefix}.fasta`: Extracted sequences from BAM file for the specified region

### Evaluation
- `{prefix}.tsv`: Evaluation metrics with QV (Quality Value) scores
- Detailed comparison between truth and query sequences with optimal matching

## Project Structure

```
haplograph/
├── Cargo.toml              # Project configuration and dependencies
├── src/
│   ├── main.rs             # Command-line interface
│   ├── asm.rs              # Assembly and graph traversal
│   ├── call.rs             # Variant calling and VCF generation
│   ├── eval.rs             # Evaluation and comparison
│   ├── extract.rs          # Sequence extraction from BAM
│   ├── graph.rs            # Graph construction
│   ├── hap.rs              # Haplotype extraction and VCF generation
│   ├── intervals.rs        # Genomic interval processing
│   ├── methyl.rs           # Methylation analysis
│   └── util.rs             # Utility functions
├── .github/
│   └── workflows/          # GitHub Actions CI/CD workflows
│       ├── ci.yml          # Continuous integration
│       ├── test-small.yml  # Quick integration tests
│       └── release.yml     # Release builds
├── wdl/                    # Workflow Definition Language files
├── scripts/                # Utility scripts
│   └── test-local.sh       # Local testing script
├── test/                   # Test data and integration tests
├── example/                # Example data and scripts
├── docker/                 # Dockerfiles for various tools
└── README.md               # This file
```

## Dependencies

### Core Dependencies
- **rust-htslib**: HTSlib bindings for BAM/VCF file handling
- **bio**: Bioinformatics algorithms and data structures
- **clap**: Command-line argument parsing
- **anyhow**: Error handling
- **serde/serde_json**: Serialization for GFA annotations

### Analysis Dependencies
- **rayon**: Parallel processing
- **ndarray**: Numerical computing
- **flate2**: Compression support
- **regex**: Pattern matching
- **csv**: CSV file processing
- **minimap2**: Sequence alignment for evaluation
- **itertools**: Iterator utilities
- **intervals**: Genomic interval operations
- **statrs**: Statistical functions
- **adjustp**: Statistical adjustments

### Development Dependencies
- **criterion**: Benchmarking
- **log/env_logger**: Logging
- **indicatif**: Progress bars

## Advanced Usage

### Custom Parameters for Different Use Cases

**High-coverage data:**
```bash
haplograph haplograph \
    --alignment-bam input.bam \
    --reference-fa reference.fa \
    --sampleid SAMPLE001 \
    --locus chr6:29943661-29943700 \
    --var-frequency-min 0.01 \
    --min-reads 5 \
    --window-size 200 \
    --primary-only
```

**Low-coverage data:**
```bash
haplograph haplograph \
    --alignment-bam input.bam \
    --reference-fa reference.fa \
    --sampleid SAMPLE001 \
    --locus chr6:29943661-29943700 \
    --var-frequency-min 0.05 \
    --min-reads 2 \
    --window-size 100
```

**Complex regions with methylation (e.g., HLA):**
```bash
haplograph haplograph \
    --alignment-bam input.bam \
    --reference-fa reference.fa \
    --sampleid SAMPLE001 \
    --locus chr6:29943661-29943700 \
    --var-frequency-min 0.02 \
    --min-reads 2 \
    --window-size 100 \
    --threshold-methyl-likelihood 0.5 \
    --file-format gfa
```

**Nanopore sequencing:**
```bash
haplograph haplograph \
    --alignment-bam input.bam \
    --reference-fa reference.fa \
    --sampleid SAMPLE001 \
    --locus chr6:29943661-29943700 \
    --detection-technology nanopore \
    --threshold-methyl-likelihood 0.5
```

## Troubleshooting

### Common Issues

1. **Memory usage**: For large regions, consider reducing window size
2. **No variants called**: Check read coverage and adjust frequency thresholds
3. **Graph assembly fails**: Ensure sufficient read overlap and adjust parameters
4. **Methylation data not found**: Ensure BAM file contains base modification tags (MM/ML)
5. **TMPDIR write errors**: Set `TMPDIR` environment variable to a writable directory (see WDL workflows)

### Performance Tips

- Use `--primary-only` for faster processing
- Adjust `--window-size` based on region complexity
- Use appropriate `--var-frequency-min` for your data type
- For interval-based analysis, use `haplointervals` instead of `haplograph`
- Consider using WDL workflows for large-scale batch processing

## Development

### Running Tests
```bash
# Unit tests
cargo test
cargo test -- --nocapture  # With output

# Local integration tests
./scripts/test-local.sh
```

### Building Documentation
```bash
cargo doc --open
```

### Code Quality
```bash
# Format code
cargo fmt --all

# Lint code
cargo clippy --all-targets --all-features -- -D warnings

# Build release
cargo build --release
```

### Continuous Integration

This project uses GitHub Actions for automated testing and building:

- **CI Workflow**: Runs on every push/PR to main/master/develop branches
  - Code formatting checks
  - Linting with clippy
  - Build verification
  - Integration tests

- **Test Workflow**: Quick integration tests for fast feedback

- **Release Workflow**: Builds release binaries for Linux, macOS, and Windows

See `.github/workflows/README.md` for more details.

## Citation

If you use Haplograph in your research, please cite:

```
Haplograph: A bioinformatics tool for haplotype analysis
Version 0.1.0
https://github.com/broadinstitute/haplograph
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## WDL Workflows

The repository includes WDL (Workflow Definition Language) files in the `wdl/` directory for running haplograph in cloud environments (e.g., Terra, Cromwell). These workflows handle:

- Batch processing of multiple samples
- Parallel execution across intervals
- Automatic resource management
- Integration with other bioinformatics tools

See the `wdl/` directory for available workflows.

## Support

For questions and support, please open an issue on the [GitHub repository](https://github.com/broadinstitute/haplograph) or contact the maintainers.
