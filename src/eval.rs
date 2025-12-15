use crate::util;
use anyhow::Result as AnyhowResult;
use bio::io::fasta::Reader as FastaReader;
use bio::io::fasta::Record;
use itertools::Itertools;
use log::info;
use minimap2::Aligner;
use rust_htslib::bam::Read as BamRead;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

pub fn find_alignment_intervals(interval_list: Vec<&str>) -> AnyhowResult<(usize, usize)> {
    //
    let mut left_bound = usize::MAX;
    let mut right_bound = usize::MIN;
    for interval in interval_list.iter() {
        let interval_ = interval.split(".").collect::<Vec<_>>()[1];
        let (chromosome, start, end) = util::split_locus(interval_.to_string());
        left_bound = left_bound.min(start);
        right_bound = right_bound.max(end);
    }
    Ok((left_bound, right_bound))
}

pub fn select_seqs(fasta: &PathBuf, haplotype_number: usize) -> AnyhowResult<Vec<Record>> {
    let fasta_records = FastaReader::from_file(fasta)?;
    // let header_list = &fasta_records.records().map(|x| x.unwrap().id().to_string()).collect::<Vec<_>>();

    // find the longest n haplotypes in the fasta file
    let mut seqs_len_vec = Vec::new();
    for record in fasta_records.records() {
        let record_ = record.unwrap();
        let record_id = record_.id().to_string();
        let region = record_id.split("\t").collect::<Vec<_>>()[0];
        let region_ = region.split(".").collect::<Vec<_>>()[0];
        let region__ = region_.split(":").collect::<Vec<_>>()[1];
        let region___ = region__.split("-").collect::<Vec<_>>();
        let left_bound = region___[0].parse::<usize>().unwrap();
        let right_bound = region___[1].parse::<usize>().unwrap();

        seqs_len_vec.push((record_.clone(), right_bound - left_bound));
    }
    //sort seqs_len_vec by the length of the sequence
    seqs_len_vec.sort_by(|a, b| b.1.cmp(&a.1));
    let seqs_len_vec = seqs_len_vec
        .iter()
        .take(haplotype_number)
        .collect::<Vec<_>>();
    let seqs_selected = seqs_len_vec.iter().map(|x| x.0.clone()).collect::<Vec<_>>();

    Ok(seqs_selected)
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
    if hits.clone().unwrap().len() != 1 {
        return Ok(-1.0);
    }
    let editdistance = hits.clone().unwrap()[0]
        .clone()
        .alignment
        .clone()
        .unwrap()
        .nm as f64;
    // let block_length = hits.clone().unwrap()[0].clone().block_len;
    let match_length = hits.clone().unwrap()[0].clone().match_len as f64;
    let alignment_length = editdistance + match_length;
    // println!("block_length == (alignment_length as i32): {}", block_length == (alignment_length as i32));
    let qv_score1 = -10_f64 * (editdistance.max(0.5) / alignment_length).log10();
    // let qv_score2 = -10 as f64* (editdistance.max(0.5) / block_length as f64).log10();
    // println!("version1: {:?}, version2: {}", qv_score1, qv_score2);

    Ok(qv_score1)
}

pub fn start(
    truth_fasta: &PathBuf,
    query_fasta: &PathBuf,
    haplotype_number: usize,
    output_prefix: &PathBuf,
) -> AnyhowResult<Vec<f64>> {
    let fasta_records = FastaReader::from_file(truth_fasta)?;
    let truth_seqs = fasta_records
        .records()
        .map(|x| x.unwrap())
        .collect::<Vec<_>>();
    let query_fasta_records = FastaReader::from_file(query_fasta)?;
    let query_seqs = query_fasta_records
        .records()
        .map(|x| x.unwrap())
        .collect::<Vec<_>>();
    let hap_num = query_seqs.len().min(haplotype_number);
    // let query_seqs_selected = select_seqs(query_fasta, hap_num)?;

    let mut optimal_score = f64::MIN;
    let mut optimal_perm = Vec::new();
    let mut optimal_qv_scores = Vec::new();

    if query_seqs.is_empty() {
        // write the evaluation results to a file
        let mut file = File::create(output_prefix)?;
        writeln!(file)?;
        return Ok(optimal_qv_scores);
    }

    let mut qv_scores_perm = Vec::new();
    for (i, t_seq) in truth_seqs.iter().enumerate() {
        for (j, q_seq) in query_seqs.iter().enumerate() {
            let qv_score = calculate_qv_score(t_seq.clone(), q_seq.clone()).unwrap();
            println!("{} {} qv_score: {}", i, j, qv_score.clone());
            qv_scores_perm.push((i, j, qv_score));
        }
    }

    for perm in qv_scores_perm.iter().permutations(hap_num) {
        let truth_seq_index = perm.iter().map(|x| x.0).collect::<Vec<_>>();
        let query_seq_index = perm.iter().map(|x| x.1).collect::<Vec<_>>();
        // remove duplicates
        if truth_seq_index.iter().unique().count() < truth_seq_index.len() {
            continue;
        }
        if query_seq_index.iter().unique().count() < query_seq_index.len() {
            continue;
        }
        let qv_scores_perm = perm.iter().map(|x| x.2).collect::<Vec<_>>();
        if qv_scores_perm.iter().sum::<f64>() > optimal_score {
            optimal_score = qv_scores_perm.iter().sum::<f64>();
            optimal_perm = perm.clone();
            optimal_qv_scores = qv_scores_perm.clone();
        }
    }

    info!("optimal_score: {}", optimal_score);
    info!("optimal_perm: {} pairs {:?}", hap_num, optimal_perm);
    info!("optimal_qv_scores: {:?}", optimal_qv_scores);

    let mut optimal_sequence_pairs = Vec::new();
    for (i, j, score) in optimal_perm.iter() {
        optimal_sequence_pairs.push((truth_seqs[*i].clone(), query_seqs[*j].clone(), *score));
    }

    // write the evaluation results to a file
    let mut file = File::create(output_prefix)?;
    for (i, j, score) in optimal_sequence_pairs.iter() {
        writeln!(file, "{} {} {}", j.id(), j.id(), score)?;
    }

    Ok(optimal_qv_scores)
}
