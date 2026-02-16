use anyhow::Result;
use std::process::Command;
use log::info;
use std::path::PathBuf;
use std::fs::File;
use std::io::{BufRead, BufReader};
use bio::io::fasta::Reader as FastaReader;
use flate2::read::GzDecoder;
use std::collections::{HashMap, HashSet};
use crate::util;
use rust_htslib::bam::Read as BamRead;
use minimap2::Aligner;
use rust_htslib::bam::record::{Cigar, CigarString};
use ndarray::Array2;
use nalgebra::DMatrix;
use std::f64;
use rayon::prelude::*;


fn cigar_to_cigarstr(cigar: &Vec<(u32, u8)>) -> CigarString {
    let op_vec: Vec<Cigar> = cigar
        .to_owned()
        .iter()
        .map(|(len, op)| match op {
            0 => Cigar::Match(*len),
            1 => Cigar::Ins(*len),
            2 => Cigar::Del(*len),
            3 => Cigar::RefSkip(*len),
            4 => Cigar::SoftClip(*len),
            5 => Cigar::HardClip(*len),
            6 => Cigar::Pad(*len),
            7 => Cigar::Equal(*len),
            8 => Cigar::Diff(*len),
            _ => panic!("Unexpected cigar operation"),
        })
        .collect();
    CigarString(op_vec)
}


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
    let records: Vec<_> = fasta_reader
        .records()
        .map(|result| result.expect("Failed to read FASTA record"))
        .collect();

    let results: Vec<_> = records
        .par_iter()
        .map(|record| {
            let seq_id = record.id().to_string();
            let seq_description = record.desc().unwrap_or("").to_string();
            let allele_id = seq_description.split(' ').next().unwrap_or("").to_string();
            let sequence = String::from_utf8_lossy(record.seq()).to_string();
            let kmer_count = count_kmer(&sequence, rollingkmer_list);
            (seq_id, allele_id, kmer_count)
        })
        .collect();

    let mut kmer_dict = HashMap::new();
    let mut seq_info = HashMap::new();
    for (seq_id, allele_id, kmer_count) in results {
        // info!("Processing sequence: {} (allele: {})", seq_id, allele_id);
        seq_info.insert(seq_id.clone(), allele_id);
        kmer_dict.insert(seq_id, kmer_count);
    }

    info!("Pangenome Extraction Completed");
    Ok((kmer_dict, seq_info))
}

pub fn count_kmer_from_read(
    pangenome_kmer_dict: &HashMap<String, HashMap<String, i32>>,
    read_path: &PathBuf,
    rollingkmer_list: &Vec<usize>,
) -> Result<(HashMap<String, i32>, HashMap<String, (String, Vec<u8>)>)> {
    
    let mut read_kmer_dict: HashMap<String, i32> = HashMap::new(); // recorded k-mer count for each read
    let mut read_seq_dict: HashMap<String, (String, Vec<u8>)> = HashMap::new(); // recorded read sequence for each read overlapping with the pangenome
    // get the set of all k-mers in the pangenome
    let pangenome_kmer_set = pangenome_kmer_dict.iter().map(|(_, kmer_dict)| kmer_dict.keys().cloned().collect::<HashSet<String>>()).flatten().collect::<HashSet<String>>();
    info!("Pangenome K-mer Set: {:?}", pangenome_kmer_set.len());
    info!("Read File: {}, {:?}", read_path.display(), read_path.display().to_string().ends_with(".bam"));

    if read_path.display().to_string().ends_with(".bam") {
        let mut bam = util::open_bam_file(&read_path.display().to_string());
        let mut read_set = HashSet::new();
        let index = 0;
        let header = bam.header().to_owned();
        for tid in 0..header.target_count() {
            let len = header.target_len(tid).unwrap_or(0) as i64;
            if len == 0 {
                continue;
            }
            bam.fetch((tid, 0, len)).expect("Failed to fetch BAM records");
            for read in bam.records() {
                let record = read.expect("Failed to read BAM record");
                let read_id = String::from_utf8_lossy(record.qname()).into_owned();
                if read_set.contains(&read_id.clone()) {
                    continue
                }
                read_set.insert(read_id.clone());
                let read_seq = String::from_utf8_lossy(&record.seq().as_bytes()).into_owned();
                let read_qual = record.qual();
                let local_kmer_dict = count_kmer(&read_seq, rollingkmer_list);
                let mut total_count = 0;
                for (kmer, count) in local_kmer_dict.iter() {
                    if !pangenome_kmer_set.contains(kmer) {
                        continue;
                    }
                    total_count += count.clone();
                    *read_kmer_dict.entry(kmer.clone()).or_insert(0) += count.clone();
                }
                if total_count > 0 {
                    read_seq_dict.insert(read_id.clone(), (read_seq.clone(), read_qual.to_vec().clone()));
                }
            }
        }
    } else if read_path.display().to_string().ends_with(".fastq") || read_path.display().to_string().ends_with(".fq") {
        return Err(anyhow::anyhow!("file format: {} UnderConstruction", read_path.display()));
    }else if read_path.display().to_string().ends_with(".fa") || read_path.display().to_string().ends_with(".fasta") {
        return Err(anyhow::anyhow!("file format: {} UnderConstruction", read_path.display()));
    } else {
        return Err(anyhow::anyhow!("Unsupported file format: {}", read_path.display()));
    };


    Ok((read_kmer_dict, read_seq_dict))

}

fn calculate_matrix(
    read_kmer_dict: &HashMap<String, i32>,
    pangenome_kmer_dict: &HashMap<String, HashMap<String, i32>>,
    sample_name: &String,
) -> Result<(Array2<f64>, Vec<String>, Vec<String>)> {

    // set row index
    let read_kmer_set = read_kmer_dict.iter().map(|(kmer, _)| kmer.clone()).collect::<HashSet<String>>();
    info!("Read K-mer Set: {:?}", read_kmer_set.len());
    let mut mrow_index = read_kmer_set.iter().cloned().collect::<Vec<String>>();
    mrow_index.sort();
    let mrow_index_dict = mrow_index.iter().enumerate().map(|(i, kmer)| (kmer.clone(), i)).collect::<HashMap<String, usize>>();
    
    //set column index
    let mut reference_index = pangenome_kmer_dict.iter().map(|(hap, _)| hap.clone()).collect::<Vec<String>>();
    reference_index.sort();
    reference_index.push(sample_name.clone());
    let mcol_index_dict = reference_index.iter().enumerate().map(|(i, hap)| (hap.clone(), i)).collect::<HashMap<String, usize>>();
    info!("Matrix Shape: {:?}, {:?}", reference_index.len(), mrow_index.len());
    // set matrix
    let mut matrix = Array2::<f64>::zeros((mrow_index.len(), reference_index.len()));
    for haplotype in reference_index.iter() {
        let kmer_dict = if haplotype == sample_name {read_kmer_dict} else {pangenome_kmer_dict.get(haplotype).unwrap()};
        for (kmer, count) in kmer_dict.iter() {
            if !mrow_index_dict.contains_key(kmer) {    
                continue;
            }
            let row_index = mrow_index_dict.get(kmer).unwrap();
            let col_index = mcol_index_dict.get(haplotype).unwrap();
            matrix[[*row_index, *col_index]] = count.clone() as f64;
        }
    }

    Ok((matrix, mrow_index, reference_index))
} 

fn calculate_best_pair(
    x: &Array2<f64>,          // features x entities
    entities: &[String],
    sample_id: &String,
) -> (String, String, f64) {
    let n = entities.len();
    let t_idx = entities.iter().position(|s| s == &sample_id.clone()).expect("test not found");

    let (rows, cols) = x.dim();
    let x_slice = x
        .as_slice()
        .expect("Array2 must be contiguous to convert to nalgebra");
    let vn = DMatrix::from_row_slice(rows, cols, x_slice);

    // normalize columns in-place
    let column_norms: Vec<f64> = vn.column_iter()
        .map(|column| column.norm())
        .collect();
    let mut vn_clone = vn.clone();
    for j in 0..cols {
        let norm = if column_norms[j] == 0.0 { 1.0 } else { column_norms[j] };
        vn_clone.column_mut(j).scale_mut(1.0 / norm);
    }

    let t_vec = vn_clone.column(t_idx).into_owned();
    let t_sim = vn_clone.transpose() * &t_vec; // (n,)
    let s = vn_clone.transpose() * &vn_clone; // (n,n)

    let (best_score, best_pair) = (0..n)
        .into_par_iter()
        .filter(|&i| i != t_idx)
        .map(|i| {
            let mut local_best_score = f64::NEG_INFINITY;
            let mut local_best_pair = (i, i);
            for j in (i + 1)..n {
                if j == t_idx {
                    continue;
                }
                let score = (t_sim[i] + t_sim[j]) / (2.0 + 2.0 * s[(i, j)]).sqrt();
                if score > local_best_score {
                    local_best_score = score;
                    local_best_pair = (i, j);
                }
            }
            (local_best_score, local_best_pair)
        })
        .reduce(
            || (f64::NEG_INFINITY, (0usize, 0usize)),
            |a, b| if a.0 >= b.0 { a } else { b },
        );

    let best_distance = 1.0 - best_score;

    (
        entities[best_pair.0].clone(),
        entities[best_pair.1].clone(),
        best_distance,
    )
}

fn realign_minimap2(
    ref_fa: &PathBuf,
    read_path: &PathBuf,
    out_bam: &PathBuf
) -> anyhow::Result<()> {

    let fastq_path = PathBuf::from(format!("{}.fastq", out_bam.display().to_string().replace(".bam", "")));
    if read_path.display().to_string().ends_with(".bam") || read_path.display().to_string().ends_with(".fasta") {
        println!("Converting BAM to FASTQ: {}", read_path.display());
        println!("Command: samtools fastq -o {} {}", read_path.display().to_string(), fastq_path.display().to_string());
        let fastq_file = std::fs::File::create(&fastq_path)?;
        Command::new("samtools")
            .arg("fastq")
            .arg(read_path.display().to_string())
            .stdout(std::process::Stdio::from(fastq_file))
            .output()?;
    } else {
        return Err(anyhow::anyhow!("Unsupporte d file format: {}", read_path.display()));
    }

    let sam_path = out_bam.with_extension("sam");
    Command::new("minimap2")
        .arg("-ax")
        .arg("asm5")
        .arg("--cs")
        .arg("--end-bonus=10")
        .arg(ref_fa.display().to_string())
        .arg(fastq_path.display().to_string())
        .arg("-o")
        .arg(sam_path.display().to_string())
        .output()?;

    Command::new("samtools")
        .arg("view")
        .arg("-bS")
        .arg(sam_path.display().to_string())
        .arg("-o")
        .arg(out_bam.display().to_string())
        .output()?;
    // sort bam
    Command::new("samtools")
        .arg("sort")
        .arg(out_bam.display().to_string())
        .arg("-o")
        .arg(format!("{}.sorted.bam", out_bam.display().to_string().replace(".bam", "")))
        .output()?;

    Command::new("samtools")
        .arg("index")
        .arg(format!("{}.sorted.bam", out_bam.display().to_string().replace(".bam", "")))
        .output()?;

    Command::new("rm")
        .arg(sam_path.display().to_string())
        .output()?;

    Command::new("rm")
        .arg(out_bam.display().to_string())
        .output()?;

    Ok(())
}
// fn realign_rust_htslib(
//     ref_fa: &PathBuf,
//     read_seq_dict: &HashMap<String, (String, Vec<u8>)>,
//     out_bam: &PathBuf,
//     data_technology: &String,
// ) -> anyhow::Result<()> {

    // // Build minimap2 index
    // let aligner = if data_technology == "hifi" || data_technology == "pacbio" {
    //     println!("Building minimap2 index for hifi/pacbio");
    //     Aligner::builder()
    //         .asm5()
    //         .with_cigar()
    //         .with_index(ref_fa, None)
    //         .expect("Unable to build minimap2 index")
    // } else if data_technology == "ont" || data_technology == "nanopore" {
    //     Aligner::builder()
    //         .map_ont()
    //         .with_cigar()
    //         .with_index(ref_fa, None)
    //         .expect("Unable to build minimap2 index")
    // }else if data_technology == "ont-r10" {
    //     Aligner::builder()
    //         .lrhq()
    //         .with_cigar()
    //         .with_index(ref_fa, None)
    //         .expect("Unable to build minimap2 index")
    // }else if data_technology == "sr" {
    //     Aligner::builder()
    //         .sr()
    //         .with_cigar()
    //         .with_index(ref_fa, None)
    //         .expect("Unable to build minimap2 index")
    // }else{
    //     Aligner::builder()
    //     .map_pb()
    //     .with_cigar()
    //     .with_index(ref_fa, None)
    //     .expect("Unable to build minimap2 index")
    // };

    // Build Header
    // Build BAM header from reference FASTA
    // let mut header = Header::new();
    // //use populate header from reference FASTA file
    // let mut ref_fa_reader = fasta::Reader::from_file(ref_fa)?;
    // for rec in ref_fa_reader.records() {
    //     let rec = rec?;
    //     let mut sq = HeaderRecord::new(b"SQ");
    //     sq.push_tag(b"SN", &rec.id());
    //     sq.push_tag(b"LN", &rec.seq().len());
    //     header.push_record(&sq);
    // }
    // let header_view = HeaderView::from_header(&header);

    // let mut writer = Writer::from_path(out_bam, &header, Format::Bam)?;

    // Align reads and write BAM
    // convert bam to fastq

    // let results = read_seq_dict.par_iter().map(|(read_id, (read_seq, qual))| {
    //     let alns = aligner
    //         .map(read_seq.as_bytes(), true, true, None, None, Some(read_id.as_bytes()))
    //         .unwrap();
    //     alns.iter()
    //         .map(|aln| {
    //             aln.to_sam_record(header_view)
    //         })
    //         .collect::<Vec<_>>()
    // }).collect::<Vec<_>>();

    // // sort results by alignment position
    // let mut final_results = results.iter().flatten().collect::<Vec<_>>();
    // final_results.sort_by_key(|r| (r.tid(), r.pos()));
    // for record in final_results.iter() {
    //     writer.write(&record).unwrap();
        
    // }

    // drop(writer);
    // index::build(out_bam, None, index::Type::Bai, 1)?;

//     Ok(())
// }

pub fn start(
    read_path: &PathBuf,
    pangenome_path: &PathBuf,
    rollingkmer_list: &Vec<usize>,
    output_prefix: &String,
    data_technology: &String,
    sample_id: &String,
) -> Result<(String, String, f64)> {
    info!("Loading Pangenome File: {}", pangenome_path.display());
    let (pangenome_kmer_dict, seq_info) = load_pangenome(pangenome_path, rollingkmer_list)?;
    info!("Counting K-mers from Read File: {}", read_path.display());
    let (read_kmer_dict, read_seq_dict) = count_kmer_from_read(&pangenome_kmer_dict, read_path, rollingkmer_list)?;
    info!("K-mer Counting Completed");
    let (matrix, mrow_index, reference_index) = calculate_matrix(&read_kmer_dict, &pangenome_kmer_dict, &sample_id)?;

    info!("Calculating Best Path");
    let best_path = calculate_best_pair(&matrix, &reference_index, sample_id);
    info!("Best Path Calculated: {:?}, {}, {}", seq_info[&best_path.0], seq_info[&best_path.1], best_path.2);
    let mut final_sequences_dict = HashMap::new();
    let fasta_reader = FastaReader::from_file(&pangenome_path.display().to_string())?;
    for record in fasta_reader.records() {
        let record = record.expect("Failed to read FASTA record");
        let seq_id = record.id().to_string();
        if seq_id == best_path.0 || seq_id == best_path.1 {
            final_sequences_dict.insert(seq_id.replace(":", "_"), String::from_utf8_lossy(record.seq()).to_string());
        }
    }
    info!("Writing Best Path to File");
    let tmp_fasta_path = PathBuf::from(format!("{}.tmp.ref.fasta", output_prefix));
    util::write_fasta(&final_sequences_dict, &tmp_fasta_path)?;
    let tmp_bam_path = PathBuf::from(format!("{}.tmp.bam", output_prefix));
    realign_minimap2(&tmp_fasta_path, &read_path, &tmp_bam_path)?;
    info!("Realigning Reads to Best Path");

    Ok(best_path)
}