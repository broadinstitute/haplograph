use anyhow::{Context, Result};
use log::info;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use bio::io::fasta::Reader as FastaReader;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use crate::util;

pub fn start(
    read_path: &PathBuf,
    pangenome_path: &PathBuf,
    rollingkmer_list: &Vec<usize>,
    output_path: &String,
) -> Result<()> {
    // Parallelize window processing - each thread gets its own BAM reader
    info!("Open Pangenome File: {}", pangenome_path.display());
    // Open FASTA file (supports both regular and gzipped files)
    let file = File::open(pangenome_path).expect("Failed to open FASTA file");
    let reader: Box<dyn BufRead> = if pangenome_path.ends_with(".gz") {
        let gz_decoder = GzDecoder::new(file);
        Box::new(BufReader::new(gz_decoder))
    } else {
        Box::new(BufReader::new(file))
    };

    let fasta_reader = FastaReader::new(reader);
    let mut kmer_dict = HashMap::new();

    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        info!("Processing sequence: {} ", seq_id);
        let sequence = String::from_utf8_lossy(record.seq()).to_string();
        let mut kmer_count = HashMap::new();
        for k in rollingkmer_list {
            for i in 0..sequence.len() - k + 1 {
                let kmer_ = &sequence[i..i+k];
                let kmer = kmer_.to_uppercase();
                let rev = util::reverse_complement(&kmer);
                *kmer_count.entry(kmer.clone()).or_insert(0) += 1;
                *kmer_count.entry(rev.clone()).or_insert(0) += 1;
            }
        }
        kmer_dict.insert(seq_id, kmer_count);
        
    }

    info!("Graph reconstruction completed");
    Ok(())
}
