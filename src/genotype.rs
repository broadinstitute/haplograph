use anyhow::{Context, Result};
use log::info;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use bio::io::fasta::Reader as FastaReader;
use flate2::read::GzDecoder;
use std::collections::{HashMap, HashSet};
use crate::util;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, Record};
use bio::io::fastq::Reader as FastqReader;
use ndarray::{Array1, Array2, Axis};
use std::f64;


pub fn count_kmer(sequence: &String, rollingkmer_list: &Vec<usize>) -> HashMap<String, i32> {
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
    kmer_count
}

pub fn load_pangenome(
    pangenome_path: &PathBuf,
    rollingkmer_list: &Vec<usize>,
) -> Result<(HashMap<String, HashMap<String, i32>>, HashMap<String, String>)> {
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
    let mut seq_info = HashMap::new();

    for result in fasta_reader.records() {
        let record = result.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        let seq_description = record.desc().unwrap_or("").to_string();
        let allele_id = seq_description.split(" ").nth(0).unwrap_or("").to_string();
        seq_info.insert(seq_id.clone(), allele_id.clone());
        info!("Processing sequence: {} (allele: {})", seq_id, allele_id);
        let sequence = String::from_utf8_lossy(record.seq()).to_string();
        let kmer_count = count_kmer(&sequence, rollingkmer_list);
        kmer_dict.insert(seq_id.clone(), kmer_count.clone());
    }

    info!("Pangenome Extraction Completed");
    Ok((kmer_dict, seq_info))
}

pub fn count_kmer_from_read(
    pangenome_kmer_dict: &HashMap<String, HashMap<String, i32>>,
    read_path: &PathBuf,
    rollingkmer_list: &Vec<usize>,
) -> Result<HashMap<String, i32>> {
    
    let mut read_kmer_dict = HashMap::new();
    // get the set of all k-mers in the pangenome
    let pangenome_kmer_set = pangenome_kmer_dict.iter().map(|(_, kmer_dict)| kmer_dict.keys().cloned().collect::<HashSet<String>>()).flatten().collect::<HashSet<String>>();
    println!("Pangenome K-mer Set: {:?}", pangenome_kmer_set);

    if read_path.ends_with(".bam") {
        let mut bam = util::open_bam_file(&read_path.display().to_string());
        for read in bam.records() {
            let record = read.expect("Failed to read BAM record");
            let read_seq = String::from_utf8_lossy(&record.seq().as_bytes()).into_owned();

            let read_kmer_count = count_kmer(&read_seq, rollingkmer_list);
            for (kmer, count) in read_kmer_count.iter() {
                if !pangenome_kmer_set.contains(kmer) {
                    continue;
                }
                if read_kmer_dict.contains_key(kmer) {
                    *read_kmer_dict.get_mut(kmer).unwrap() += count;
                } else {
                    read_kmer_dict.insert(kmer.clone(), count.clone());
                }
            }
        }
    } else if read_path.ends_with(".fastq") || read_path.ends_with(".fq") {
        let fastq_reader = FastqReader::from_file(&read_path.display().to_string())?;
        for record in fastq_reader.records() {
            let record = record.expect("Failed to read FASTQ record");
            
            let read_seq = String::from_utf8_lossy(record.seq()).to_string();
            let read_kmer_count = count_kmer(&read_seq, rollingkmer_list);
            for (kmer, count) in read_kmer_count.iter() {
                if !pangenome_kmer_set.contains(kmer) {
                    continue;
                }
                if read_kmer_dict.contains_key(kmer) {
                    *read_kmer_dict.get_mut(kmer).unwrap() += count;
                } else {
                    read_kmer_dict.insert(kmer.clone(), count.clone());
                }
            }
        }

    }else if read_path.ends_with(".fa") || read_path.ends_with(".fasta") {
        let fasta_reader = FastaReader::from_file(&read_path.display().to_string())?;
        for record in fasta_reader.records() {
            let record = record.expect("Failed to read FASTA record");
            
            let read_seq = String::from_utf8_lossy(record.seq()).to_string();
            let read_kmer_count = count_kmer(&read_seq, rollingkmer_list);
            for (kmer, count) in read_kmer_count.iter() {
                if !pangenome_kmer_set.contains(kmer) {
                    continue;
                }
                if read_kmer_dict.contains_key(kmer) {
                    *read_kmer_dict.get_mut(kmer).unwrap() += count;
                } else {
                    read_kmer_dict.insert(kmer.clone(), count.clone());
                }
            }
        }

    }else{
        return Err(anyhow::anyhow!("Unsupported file format: {}", read_path.display()));
    }
    Ok(read_kmer_dict)
}



fn calculate_best_path(
    read_kmer_dict: &HashMap<String, i32>,
    pangenome_kmer_dict: &HashMap<String, HashMap<String, i32>>,
    seq_info: &HashMap<String, String>,
) -> Result<(String, String, f64)> {
    // set row index
    let read_kmer_set = read_kmer_dict.iter().map(|(kmer, _)| kmer.clone()).collect::<HashSet<String>>();
    let pangenome_kmer_set = pangenome_kmer_dict.iter().map(|(_, kmer_dict)| kmer_dict.keys().cloned().collect::<HashSet<String>>()).flatten().collect::<HashSet<String>>();
    let mut mrow_index = read_kmer_set.union(&pangenome_kmer_set).cloned().collect::<Vec<String>>();
    let mrow_index_dict = mrow_index.iter().enumerate().map(|(i, kmer)| (kmer.clone(), i)).collect::<HashMap<String, usize>>();
    mrow_index.sort();
    

    //set column index
    let mut reference_index = pangenome_kmer_dict.iter().map(|(hap, _)| hap.clone()).collect::<Vec<String>>();
    reference_index.sort();
    reference_index.push("test".to_string());
    let mcol_index_dict = reference_index.iter().enumerate().map(|(i, hap)| (hap.clone(), i)).collect::<HashMap<String, usize>>();

    // set matrix
    let mut matrix = Array2::<f64>::zeros((mrow_index.len(), reference_index.len()));
    for (haplotype, kmer_dict) in pangenome_kmer_dict.iter() {
        for (kmer, count) in kmer_dict.iter() {
            let row_index = mrow_index_dict.get(kmer).unwrap();
            let col_index = mcol_index_dict.get(haplotype).unwrap();
            matrix[[*row_index, *col_index]] = *count as f64;
        }
    }
    // get shape of matrix
    let n = reference_index.len();
    let t_idx = reference_index
        .iter()
        .position(|s| s == &"test".to_string())
        .expect("test not found");

    // normalize matrix by column
    let mut norms = Array1::<f64>::zeros(n);
    for (j, col) in matrix.axis_iter(Axis(1)).enumerate() {
        let norm = col.dot(&col).sqrt();
        norms[j] = if norm == 0.0 { 1.0 } else { norm };
    }
    // normalize matrix by column
    for (j, mut col) in matrix.axis_iter_mut(Axis(1)).enumerate() {
        col /= norms[j];
    }

    // t_vec, t_sim, and S
    let t_vec = matrix.column(t_idx).to_owned();
    let t_sim = matrix.t().dot(&t_vec);      // shape (n,)
    let s = matrix.t().dot(&matrix);   

    // brute-force all pairs (i<j, i/j != test)
    let mut best_score = f64::NEG_INFINITY;
    let mut best_pair = (0usize, 0usize);

    for i in 0..n {
        if i == t_idx { continue; }
        for j in (i + 1)..n {
            if j == t_idx { continue; }
            let den = (2.0 + 2.0 * s[[i, j]]).sqrt();
            let num = t_sim[i] + t_sim[j];
            let score = num / den;
            if score > best_score {
                best_score = score;
                best_pair = (i, j);
            }
        }
    } 


    let best_distance = 1.0 - best_score;
    let best_path_1 = reference_index[best_pair.0].clone();
    let best_path_2 = reference_index[best_pair.1].clone();


    Ok((
        best_path_1,
        best_path_2,
        best_distance,
    ))
} 


pub fn start(
    read_path: &PathBuf,
    pangenome_path: &PathBuf,
    rollingkmer_list: &Vec<usize>,
    output_prefix: &String,
) -> Result<(String, String, f64)> {
    info!("Loading Pangenome File: {}", pangenome_path.display());
    let (pangenome_kmer_dict, seq_info) = load_pangenome(pangenome_path, rollingkmer_list)?;
    info!("Counting K-mers from Read File: {}", read_path.display());
    let read_kmer_dict = count_kmer_from_read(&pangenome_kmer_dict, read_path, rollingkmer_list)?;
    info!("K-mer Counting Completed");

    info!("Calculating Best Path");
    let best_path = calculate_best_path(&read_kmer_dict, &pangenome_kmer_dict, &seq_info)?;
    let mut final_sequences_dict = HashMap::new();
    let fasta_reader = FastaReader::from_file(&pangenome_path.display().to_string())?;
    for record in fasta_reader.records() {
        let record = record.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        if seq_id == best_path.0 || seq_id == best_path.1 {
            final_sequences_dict.insert(seq_id, String::from_utf8_lossy(record.seq()).to_string());
        }
    }
    util::write_graph_path_fasta(&final_sequences_dict, &PathBuf::from(format!("tmp.{}.fasta", output_prefix)))?;
    Ok(best_path)
}