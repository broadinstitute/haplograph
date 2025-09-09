use std::collections::{HashMap, HashSet};
use url::Url;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use log::{info, warn};
use anyhow::{Result as AnyhowResult};
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use bio::io::fasta::Reader as FastaReader;
use bio::io::fastq;

pub fn select_seqs(fasta: &FastaReader, haplotype_number: usize) -> AnyhowResult<(Vec<fastq::Record>)> {
    // find the longest n haplotypes in the fasta file
    let mut seqs_len_vec = Vec::new();
    for record in fasta.records(){
        seqs_len_vec.push((record.id().to_string(), record.seq().len()));
    }
    //sort seqs_len_vec by the length of the sequence
    seqs_len_vec.sort_by(|a, b| b.1.cmp(&a.1));
    let seqs_len_vec = seqs_len_vec.iter().take(haplotype_number).collect::<Vec<_>>();
    let seqs_set = seqs_len_vec.iter().map(|x| x.0.to_string()).collect::<HashSet<_>>();
    let mut seqs_selected = Vec::new();
    for record in fasta.records(){
        if seqs_set.contains(&record.id().to_string()){
            seqs_selected.push(record.clone());
        }
    }
    Ok((seqs_selected))
}

pub fn start(truth_fasta: &PathBuf, query_fasta: &Vec<fastq::Record>, output_prefix: &PathBuf) -> AnyhowResult<()> {
    let truth_fasta = FastaReader::from_path(truth_fasta)?;
    let query_fasta = FastaReader::from_path(query_fasta)?;
    let truth_seqs = truth_fasta.records().collect::<Vec<_>>();
    let query_seqs = query_fasta.records().collect::<Vec<_>>();
    let truth_seqs_len = truth_seqs.len();
    let query_seqs_len = query_seqs.len();
    let mut truth_seqs_set = HashSet::new();
    let mut query_seqs_set = HashSet::new();
    for seq in truth_seqs {
    Ok(())
}