# Haplograph

A bioinformatics tool for haplotype analysis from BAM files.

## Features

- Extract aligned reads from BAM files within specified genomic regions
- Support for haplotype-specific filtering
- Command-line interface with comprehensive argument parsing
- Library API for programmatic use
- Comprehensive error handling and logging

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/haplograph.git
cd haplograph

# Build the project
cargo build --release

# Install globally (optional)
cargo install --path .
```

## Usage

### Command Line Interface

```bash
# Basic usage
haplograph -b input.bam -c chr1 -s 1000 -e 2000 -o output/

# With haplotype filtering
haplograph -b input.bam -c chr1 -s 1000 -e 2000 -h 1 -o output/

# Verbose output
haplograph -b input.bam -c chr1 -s 1000 -e 2000 -v -o output/
```

### Programmatic Usage

```rust
use haplograph::{extract_aligned_bam_reads, validate_file_path};
use rust_htslib::bam::IndexedReader;

fn main() -> anyhow::Result<()> {
    let mut bam = IndexedReader::from_path("input.bam")?;
    
    let reads = extract_aligned_bam_reads(
        "analysis",
        &mut bam,
        "chr1",
        &1000,
        &2000,
        "extracted",
        None,
    )?;
    
    println!("Extracted {} reads", reads.len());
    Ok(())
}
```

## Project Structure

```
haplograph/
├── Cargo.toml          # Project configuration and dependencies
├── src/
│   ├── main.rs         # Command-line interface entry point
│   ├── lib.rs          # Library entry point and public API
│   ├── build.rs        # Core BAM processing functionality
│   └── utils.rs        # Utility functions
├── tests/
│   └── integration_test.rs  # Integration tests
├── examples/
│   └── basic_usage.rs  # Example usage
└── README.md           # This file
```

## Dependencies

- **bio**: Bioinformatics algorithms and data structures
- **rust-htslib**: HTSlib bindings for BAM/CRAM file handling
- **clap**: Command-line argument parsing
- **anyhow**: Error handling
- **log/env_logger**: Logging
- **serde**: Serialization/deserialization
- **rayon**: Parallel processing
- **ndarray**: Numerical computing

## Development

### Running Tests

```bash
# Run unit tests
cargo test

# Run integration tests
cargo test --test integration_test

# Run with output
cargo test -- --nocapture
```

### Building Documentation

```bash
# Generate documentation
cargo doc --open

# Generate documentation for all dependencies
cargo doc --document-private-items --open
```

### Code Formatting and Linting

```bash
# Format code
cargo fmt

# Run clippy linter
cargo clippy

# Run clippy with all warnings
cargo clippy -- -W clippy::all
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## Citation

If you use this tool in your research, please cite:

```
Haplograph: A bioinformatics tool for haplotype analysis
Version 0.1.0
https://github.com/yourusername/haplograph
```