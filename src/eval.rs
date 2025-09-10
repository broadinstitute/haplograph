use std::collections::{HashMap, HashSet};
use url::Url;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use log::{info, warn};
use anyhow::{Result as AnyhowResult};
use bio::alignment::pairwise::*;
use bio::alignment::AlignmentOperation;
use bio::io::fasta::Reader as FastaReader;
use bio::io::fastq;
use bio::io::fasta::Record;
use std::path::PathBuf;
use minimap2::{Aligner, Preset};
use itertools::Itertools;


pub fn select_seqs(fasta: &PathBuf, haplotype_number: usize) -> AnyhowResult<(Vec<Record>)> {
    let fasta_records = FastaReader::from_file(fasta)?;
    // find the longest n haplotypes in the fasta file
    let mut seqs_len_vec = Vec::new();
    for record in fasta_records.records(){
        let record_ = record.unwrap();
        seqs_len_vec.push((record_.clone(), record_.seq().len()));
    }
    //sort seqs_len_vec by the length of the sequence
    seqs_len_vec.sort_by(|a, b| b.1.cmp(&a.1));
    let seqs_len_vec = seqs_len_vec.iter().take(haplotype_number).collect::<Vec<_>>();
    let seqs_selected = seqs_len_vec.iter().map(|x| x.0.clone()).collect::<Vec<_>>();

    Ok((seqs_selected))
}

pub fn calculate_qv_score(truth_seqs: Record, query_seqs: Record) -> AnyhowResult<f64> {
    let truth_seq = truth_seqs.seq();
    let query_seq = query_seqs.seq();
    let aligner = Aligner::builder()
                                            .asm5()
                                            .with_cigar()
                                            .with_seq(truth_seq)
                                            .expect("Unable to build index");
    let hits = aligner.map(query_seq, true, true, None, None, Some(b"Query Name"));
    assert_eq!(hits.clone().unwrap().len(), 1);
    let editdistance = hits.clone().unwrap()[0].clone().alignment.clone().unwrap().nm as f64;
    // let block_length = hits.clone().unwrap()[0].clone().block_len;
    let match_length = hits.clone().unwrap()[0].clone().match_len as f64;
    let alignment_length = editdistance + match_length;
    // println!("block_length == (alignment_length as i32): {}", block_length == (alignment_length as i32));
    let qv_score1 = -10 as f64* (editdistance.max(0.5) / alignment_length).log10();
    // let qv_score2 = -10 as f64* (editdistance.max(0.5) / block_length as f64).log10();
    // println!("version1: {:?}, version2: {}", qv_score1, qv_score2);

    Ok(qv_score1)
}

pub fn start(truth_fasta: &PathBuf, query_fasta: &PathBuf, haplotype_number: usize, output_prefix: &PathBuf) -> AnyhowResult<(Vec<f64>)> {
    let truth_seqs = select_seqs(truth_fasta, haplotype_number)?;
    let query_seqs = select_seqs(query_fasta, haplotype_number)?;
    println!("truth_seqs: {}", truth_seqs.len());
    println!("query_seqs: {}", query_seqs.len());

    let v = (0..haplotype_number).collect::<Vec<_>>();
    let mut optimal_score = f64::MIN;
    let mut optimal_perm = Vec::new();
    let mut optimal_qv_scores = Vec::new();
    for perm in v.iter().permutations(v.len()).unique() {
        let mut qv_scores_perm = Vec::new();
        for (i, j) in perm.iter().enumerate() {
            let qv_score = calculate_qv_score(truth_seqs[i].clone(), query_seqs[**j].clone()).unwrap();
            println!("{} {} qv_score: {}", i, j, qv_score.clone());
            qv_scores_perm.push(qv_score.clone());
        }
        if qv_scores_perm.iter().sum::<f64>() > optimal_score {
            optimal_score = qv_scores_perm.iter().sum::<f64>();
            optimal_perm = perm.clone();
            optimal_qv_scores = qv_scores_perm.clone();
        }
    }
    println!("optimal_score: {}", optimal_score);
    println!("optimal_perm: {:?}, {:?}", v, optimal_perm);
    println!("optimal_qv_scores: {:?}", optimal_qv_scores);


    Ok((optimal_qv_scores))
}