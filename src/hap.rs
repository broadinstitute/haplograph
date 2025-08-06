use log::{info, warn};
use anyhow::{Result, Context};
use rust_htslib::bam::IndexedReader;
use std::path::{PathBuf};
use url::Url;
use bio::io::fastq;
use std::collections::HashMap;
use std::io::Write;
use std::fs::File;
use rust_htslib::faidx::Reader;
use crate::Cli;
use crate::util;


fn process_bam_file(cli: &Cli, start: usize, end: usize) -> Vec<fastq::Record> {
    
    // Open BAM file
    if cli.alignment_bam.to_string().starts_with("gs://") && std::env::var("GCS_OAUTH_TOKEN").is_err() {
        util::gcs_authorize_data_access();
    }
    let mut bam = if cli.alignment_bam.starts_with("gs://") {
        let url = Url::parse(&cli.alignment_bam).unwrap();
        IndexedReader::from_url(&url).unwrap()
    } else {
        IndexedReader::from_path(&cli.alignment_bam).unwrap()
    };
    
    info!("Successfully opened BAM file");
    
    // Extract reads from the specified region
    let reads = util::extract_reads_at_single_locus(
        &mut bam,
        &cli.chromosome,
        start as u64,
        end as u64,
        Some(cli.quality), // min_mapping_quality
    ).unwrap();
    
     info!("Extracted {} reads from region", reads.len());

    reads
}

fn process_fasta_file(cli: &Cli, start:u64, end: u64) -> Vec<fastq::Record> {
    
    // Open BAM file
    if cli.reference_fa.to_string().starts_with("gs://") && std::env::var("GCS_OAUTH_TOKEN").is_err() {
        util::gcs_authorize_data_access();
    }
    let mut fasta = if cli.reference_fa.starts_with("gs://") {
        let url = Url::parse(&cli.reference_fa).unwrap();
        Reader::from_url(&url).unwrap()
    } else {
        Reader::from_path(&cli.reference_fa).unwrap()
    };
    
    info!("Successfully opened FASTA file");
    
    // Extract reads from the specified region
    util::extract_fasta_seqs(
        &cli.sampleid,
        &mut fasta,
        &cli.chromosome,
        &start,
        &end,
        &"".to_string()
        ).unwrap()
}


fn write_fasta_output(reads: &HashMap<String, Vec<String>>, output_path: &PathBuf) -> Result<()> {
    
    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create output file: {}", output_path.display()))?;
    
    let mut index:i32 = 1;

    for (sequence, read_names) in reads {
        // Use the first read name as the sequence ID, or create a hash-based ID
        let sequence_id = if !read_names.is_empty() {
            read_names.join(",")
        } else {
            "unknown".to_string()   
        };
        
        // Write FASTA format: >sequence_id\nsequence\n
        writeln!(file, ">Haplotype{}\tLength{}\tSupportReads {}\t{}", index, sequence.len(), read_names.len(), sequence_id)?;
        writeln!(file, "{}", sequence)?;
        index += 1;
    }
    
    Ok(())
}
fn extract_haplotypes (cli: &Cli, start: usize, end: usize) -> Result<HashMap<String, Vec<String>>> {
    // Validate inputs
    let reads = if start >= end {
        anyhow::bail!("Start position must be less than end position");
    } else {
        // Continue with processing
        info!("Input validation passed");
        // Process the BAM file
        process_bam_file(cli, start, end)
    };
    
    let reference = if cli.begin >= cli.end {
        anyhow::bail!("Start position must be less than end position");
    } else {
        // Continue with processing
        info!("reference validation passed");
        // Process the BAM file
        process_fasta_file(cli, start as u64, end as u64)
    };
    
    // Read the reads records (name and sequence) into a vector.
    let mut read_dictionary: HashMap<String, Vec<String>> = HashMap::new();
    for read in reads{
        let contig_name = read.id().to_string();
        let contig = String::from_utf8_lossy(read.seq()).to_string();
        read_dictionary.entry(contig).or_default().push(contig_name);
    }
    for read in reference {
        let contig_name = read.id().to_string();
        let contig = String::from_utf8_lossy(read.seq()).to_string();
        read_dictionary.entry(contig).or_default().push(contig_name);
    }
    let mut final_hap: HashMap<String, Vec<String>> = HashMap::new();
    for (hap, vec) in read_dictionary {
        if vec.len() >= cli.min_reads as usize {
            final_hap.insert(hap, vec);
        }
    }
    info!("Extracted {} haplotypes from region", final_hap.len());
    Ok(final_hap)
}

pub fn start(cli: &Cli, start: usize, end: usize) -> Result<()> {
    // check if the region is too large
    if cli.end - cli.begin > cli.limited_size {
        warn!("Region is too large, please consider using a smaller region");
    }

    let final_hap = extract_haplotypes(cli, cli.begin as usize, cli.end as usize)?;
    
    let final_hap_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.fasta", cli.output, cli.sampleid, cli.chromosome, cli.begin, cli.end));
    write_fasta_output(&final_hap, &final_hap_output)?;

    info!("Haplotype reconstruction completed");
    Ok(())
}