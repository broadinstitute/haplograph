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

pub fn process_bam_file(cli: &Cli, start: usize, end: usize) -> Vec<fastq::Record> {
    
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

pub fn process_fasta_file(cli: &Cli, start: u64, end: u64) -> Vec<fastq::Record> {
    use bio::io::fasta::Reader as FastaReader;
    
    info!("Opening FASTA file: {}", cli.reference_fa);
    
    // Open FASTA file
    let mut fasta_reader = FastaReader::from_file(&cli.reference_fa)
        .expect("Failed to open FASTA file");
    
    let mut reference_seqs = Vec::new();
    
    // Read all sequences from FASTA
    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        let sequence = String::from_utf8_lossy(record.seq()).to_string();
        
        info!("Found sequence: {} (length: {})", seq_id, sequence.len());
        
        // Check if this sequence matches the target chromosome
        if seq_id == cli.chromosome {
            info!("Extracting region {}:{}-{} from chromosome {}", 
                  cli.chromosome, start, end, seq_id);
            
            // Extract the specified region
            if start < sequence.len() as u64 && end <= sequence.len() as u64 {
                let region_seq = sequence[start as usize..end as usize].to_string();
                let record_id = format!("{}:{}-{}|reference|{}", 
                                       cli.chromosome, start, end, cli.sampleid);
                
                let fastq_record = fastq::Record::with_attrs(
                    &record_id,
                    None,
                    region_seq.as_bytes(),
                    vec![30; region_seq.len()].as_slice() // Default quality score
                );
                
                reference_seqs.push(fastq_record);
                info!("Extracted reference sequence: {} bp", region_seq.len());
            } else {
                warn!("Region coordinates out of bounds for sequence length {}", sequence.len());
            }
        }
    }
    
    if reference_seqs.is_empty() {
        warn!("No matching chromosome '{}' found in FASTA file", cli.chromosome);
    }
    
    info!("Extracted {} reference sequences", reference_seqs.len());
    reference_seqs
}


pub fn write_fasta_output(reads: &HashMap<String, Vec<String>>, output_path: &PathBuf) -> Result<()> {

    
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
        // writeln!(file, ">{}{}\t{}{}\t{}{}\t{}", "Haplotype", index, "Length", sequence.len(), "SupportReads ", read_names.len(), sequence_id)?;
        writeln!(file, ">Haplotype{}\tLength{}\tSupportReads {}\t{}", index, sequence.len(), read_names.len(), sequence_id)?;
        writeln!(file, "{}", sequence)?;
        index += 1;
    }
    
    Ok(())
}

pub fn extract_haplotypes (cli: &Cli, start: usize, end: usize) -> Result<HashMap<String, Vec<String>>> {
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
    info!("Extracted {} haplotypes from region, supported reads {:?}", final_hap.len(), final_hap.values().map(|v| v.len()).sum::<usize>());
    Ok(final_hap)
}

pub fn start(cli: &Cli, start: usize, end: usize) -> Result<()> {
    // check if the region is too large
    if cli.end - cli.begin > cli.limited_size {
        warn!("Region is too large, please consider using a smaller region");
    }

    let final_hap = extract_haplotypes(cli, start as usize, end as usize)?;
    
    let final_hap_output = PathBuf::from(format!("{}/{}_{}_{}_{}_haplograph.fasta", cli.output, cli.sampleid, cli.chromosome, start, end));
    write_fasta_output(&final_hap, &final_hap_output)?;

    info!("Haplotype reconstruction completed");
    Ok(())
}