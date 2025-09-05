use log::{info, warn};
use anyhow::{Result as AnyhowResult, Context};
use bio::io::fastq;
use rust_htslib::faidx::Reader;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, record::Aux};
use std::path::{PathBuf};
use url::Url;
use std::collections::{HashMap, BTreeMap, HashSet};
use std::io::Write;
use std::fs::File;
use crate::Cli;
use crate::util;
use crate::intervals;
use std::error::Error;
use indicatif::{ProgressBar, ProgressStyle};

pub fn start(cli: &Cli, bam: &mut IndexedReader, chromosome: &str, windows: &Vec<(usize, usize)>) -> AnyhowResult<()> {

    let (chromosome, start, end) = util::split_locus(cli.locus.clone());
    for (i, window) in windows.iter().enumerate() {
        let (start, end) = window;
        let haplotype_info = intervals::start(bam, &cli.reference_fa, &chromosome, *start, *end, &cli.sampleid, cli.min_reads as usize, cli.frequency_min, cli.primary_only, true).unwrap();
    }

    info!("Haplotype reconstruction completed");
    Ok(())
}