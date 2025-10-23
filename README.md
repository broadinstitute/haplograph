# Haplograph

A  bioinformatics tool for haplotype/meplotype analysis and variant calling from BAM files using graph-based assembly approaches.

## Overview

Haplograph is designed for accurate haplotype or meplotype (methylation aware) reconstruction and variant calling in genomic regions, particularly useful for complex region such as HLA/MHC analysis, and other highly polymorphic regions. The tool uses graph-based assembly to reconstruct haplotypes (methylation aware) and call variants with high accuracy.

## Features

- **Graph-based Haplotype Assembly**: Constructs sequence graphs from aligned reads
- **Variant Calling**: Generates VCF files with phased variants
- **Multiple Output Formats**: Supports both GFA (Graph Fragment Assembly), VCF and FASTA formats
- **Germline Haplotype Filtering**: Focus on major haplotypes for cleaner results
- **Read-based Phasing**: Uses read overlap information for accurate phasing
- **Comprehensive Evaluation**: Compare results against truth sets
- **Flexible Parameters**: Configurable window sizes, read support thresholds, and filtering options

## Installation

### Prerequisites

- Rust 1.70+ (install from [rustup.rs](https://rustup.rs/))
- HTSlib development libraries
- Standard bioinformatics tools (samtools, bcftools) for file processing

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

Haplograph provides four main commands for different stages of haplotype analysis:

### 1. Haplotype Extraction

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
- `--frequency-min`: Minimum allele frequency (default: 0.01)
- `--min-reads`: Minimum supporting reads (default: 2)
- `--window-size`: Window size for analysis (default: 1000)
- `--primary-only`: Use only primary alignments
- `--default-file-format`: Output format (fasta or gfa, default: gfa)

### 2. Haplotype Assembly

Assemble haplotypes from GFA files:

```bash
# Basic assembly
haplograph assemble \
    --graph-gfa output/HLA_A.gfa \
    --output-prefix output/HLA_A_asm

# Germline-only assembly with specific haplotype count
haplograph assemble \
    --graph-gfa output/HLA_A.gfa \
    --major-haplotype-only \
    --number-of-haplotypes 2 \
    --output-prefix output/HLA_A_asm \
    --verbose
```

**Parameters:**
- `--graph-gfa`: Input GFA file
- `--locus`: Genomic region
- `--major-haplotype-only`: Focus on major haplotypes only
- `--number-of-haplotypes`: Number of haplotypes to extract (default: 2)

### 3. Variant Calling

Call variants from assembled haplotypes:

```bash
# Basic variant calling
haplograph call \
    --gfa-file output/HLA_A_asm.gfa \
    --sampleid SAMPLE001 \
    --reference-fa reference.fa \
    --output-prefix output/HLA_A_variants \
    --phase-variants

# With phasing
haplograph call \
    --gfa-file output/HLA_A_asm.gfa \
    --sampleid SAMPLE001 \
    --reference-fa reference.fa \
    --phase-variants \
    --maximum-haplotypes 2 \
    --output-prefix output/HLA_A_variants \
    --verbose
```

**Parameters:**
- `--gfa-file`: Input GFA file from assembly
- `--sampleid`: Sample identifier
- `--reference-fa`: Reference FASTA file
- `--phase-variants`: Enable variant phasing
- `--maximum-haplotypes`: Maximum number of haplotypes (default: 2)

### 4. Evaluation

Evaluate haplotype accuracy against truth sets:

```bash
haplograph evaluate \
    --truth-fasta truth_haplotypes.fasta \
    --query-fasta output/HLA_A_asm.fasta \
    --seq-number 2 \
    --output-prefix output/evaluation \
    --verbose
```

**Parameters:**
- `--truth-fasta`: Truth haplotype sequences
- `--query-fasta`: Query haplotype sequences
- `--seq-number`: Number of sequences to compare (default: 2)

## Complete Workflow Example

```bash
# 1. Extract haplotypes and build graph
haplograph haplograph \
    --alignment-bam sample.bam \
    --reference-fa hg38.fa \
    --sampleid HG002 \
    --locus chr6:29943661-29943700 \
    --output-prefix output/HLA_A

# 2. Assemble haplotypes
haplograph assemble \
    --graph-gfa output/HLA_A.gfa \
    --major-haplotype-only \
    --number-of-haplotypes 2 \
    --output-prefix output/HLA_A_asm

# 3. Call variants
haplograph call \
    --gfa-file output/HLA_A_asm.gfa \
    --sampleid HG002 \
    --reference-fa hg38.fa \
    --phase-variants \
    --output-prefix output/HLA_A_variants

# 4. Evaluate results (if truth available)
haplograph evaluate \
    --truth-fasta truth_HLA_A.fasta \
    --query-fasta output/HLA_A_asm.fasta \
    --output-prefix output/evaluation
```

## Output Files

### Haplotype Extraction
- `{prefix}.gfa`: Sequence graph in GFA format
- `{prefix}.fasta`: Haplotype sequences (if fasta format selected)

### Assembly
- `{prefix}.fasta`: Assembled haplotype sequences
- Contains multiple haplotype sequences with support information

### Variant Calling
- `{prefix}.vcf.gz`: Compressed VCF file with called variants
- `{prefix}.vcf.gz.tbi`: Tabix index for the VCF file
- Includes phased variants with quality scores

### Evaluation
- `{prefix}.tsv`: Evaluation metrics and accuracy scores
- Detailed comparison between truth and query sequences

## Project Structure

```
haplograph/
├── Cargo.toml              # Project configuration and dependencies
├── src/
│   ├── main.rs             # Command-line interface
│   ├── asm.rs              # Assembly and graph traversal
│   ├── call.rs             # Variant calling and VCF generation
│   ├── eval.rs             # Evaluation and comparison
│   ├── graph.rs            # Graph construction
│   ├── hap.rs              # Haplotype extraction
│   ├── intervals.rs        # Genomic interval processing
│   └── util.rs             # Utility functions
├── wdl/                    # Workflow Definition Language files
├── example/                # Example data and scripts
├── output/                 # Output directory (created during analysis)
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

### Development Dependencies
- **criterion**: Benchmarking
- **log/env_logger**: Logging
- **indicatif**: Progress bars

## Advanced Usage

### Custom Parameters for Different Use Cases

**High-coverage data:**
```bash
haplograph haplotype \
    --frequency-min 0.01 \
    --min-reads 5 \
    --window-size 200 \
    --primary-only
```

**Low-coverage data:**
```bash
haplograph haplotype \
    --frequency-min 0.05 \
    --min-reads 2 \
    --window-size 100
```

**Complex regions (e.g., HLA):**
```bash
haplograph haplotype \
    --frequency-min 0.02 \
    --min-reads 2 \
    --window-size 100 \
    --default-file-format gfa
```

## Troubleshooting

### Common Issues

1. **Memory usage**: For large regions, consider reducing window size
2. **No variants called**: Check read coverage and adjust frequency thresholds
3. **Graph assembly fails**: Ensure sufficient read overlap and adjust parameters

### Performance Tips

- Use `--primary-only` for faster processing
- Adjust `--window-size` based on region complexity
- Use appropriate `--frequency-min` for your data type

## Development

### Running Tests
```bash
cargo test
cargo test -- --nocapture  # With output
```

### Building Documentation
```bash
cargo doc --open
```

### Code Quality
```bash
cargo fmt
cargo clippy
```

## Citation

If you use Haplograph in your research, please cite:

```
Haplograph: A bioinformatics tool for haplotype analysis
Version 0.1.0
https://github.com/broadinstitute/haplograph
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Support

For questions and support, please open an issue on the GitHub repository or contact the maintainers.
