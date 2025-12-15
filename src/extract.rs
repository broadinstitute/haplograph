use crate::intervals;
use anyhow::{Context, Result as AnyhowResult};
use rust_htslib::bam::IndexedReader;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

pub fn start(
    bam: &mut IndexedReader,
    chromosome: &str,
    start: usize,
    end: usize,
    primary_only: bool,
    output_path: String,
) -> AnyhowResult<()> {
    let (reads, read_coordinates, read_sequence_dictionary, bam_records) =
        intervals::extract_haplotypes_coordinates_from_bam(
            bam,
            chromosome,
            start as u64,
            end as u64,
            primary_only,
        )
        .unwrap();
    let outputfile = PathBuf::from(format!("{}.fasta", output_path));
    let mut file = File::create(outputfile)
        .with_context(|| format!("Failed to create output file: {}.fasta", output_path))?;

    for r in reads.iter() {
        let r_name = String::from_utf8_lossy(r.id().as_bytes()).to_string();
        let r_seq = String::from_utf8_lossy(r.seq()).to_string();
        writeln!(file, ">{}", r_name)?;
        // write the sequence in fasta format
        let seq_len = r_seq.len();
        let chars_per_line = 60;
        let full_lines = seq_len / chars_per_line;
        for i in 0..full_lines {
            let start = i * chars_per_line;
            let end = start + chars_per_line;
            writeln!(file, "{}", &r_seq[start..end])?;
        }
        if seq_len % chars_per_line != 0 {
            writeln!(file, "{}", &r_seq[full_lines * chars_per_line..])?;
        }
    }
    Ok(())
}
