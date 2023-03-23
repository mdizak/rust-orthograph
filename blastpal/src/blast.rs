use crate::models::{Blast, HmmSearch};
use biotools::settings::Settings;
use lazy_static::lazy_static;
use log::info;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::io::{self, BufRead};
use std::path::Path;
use std::process::Command;

lazy_static! {
    pub static ref CONFIG: Settings = Settings::new();
}

pub fn run(search: &HmmSearch) -> Vec<Blast> {
    // Generate tmp file
    let tmpfile = generate_tmpfile(&search);

    // Get filename
    let ali_prefix = format!("{{{0}-{1}}}", (search.ali_start - 1), (search.ali_end - 1));
    let outfile = format!(
        "{}/blast/{}-{}-{}-{}.blast",
        CONFIG.report.output_dir, search.gene_id, search.target, ali_prefix, CONFIG.report.set_name
    );

    // Delete file, if exists
    if Path::new(&outfile).exists() {
        fs::remove_file(&outfile).expect("Unable to delete previous blast results file");
    }

    // Get blastdb
    let blastdb = if CONFIG.report.blastdb == "" {
        format!(
            "{}/{}/blast/{}",
            CONFIG.report.sets_dir, CONFIG.report.set_name, CONFIG.report.set_name
        )
    } else {
        format!("{}", CONFIG.report.blastdb)
    };

    // Run blastp command
    let output = Command::new(&CONFIG.programs.blast)
        .args([
            "-outfmt",
            "7 qseqid sseqid evalue bitscore qstart qend",
            "-evalue",
            &CONFIG.search.blast_evalue_threshold.to_string(),
            "-threshold",
            &CONFIG.search.blast_threshold.to_string(),
            "-max_target_seqs",
            &CONFIG.search.max_blast_searches.to_string(),
            "-num_threads",
            &CONFIG.search.num_threads.to_string(),
            "-db",
            &blastdb,
            "-query",
            &tmpfile,
            "-out",
            &outfile,
        ])
        .output()
        .expect("Unable to run blastp program");

    // Log
    info!(
        "Completed blast of hmm search id# {} with {} at file: {}",
        search.id, output.status, outfile
    );
    info!("blastp stdout: {}", String::from_utf8_lossy(&output.stdout));
    info!("blastp stderr: {}", String::from_utf8_lossy(&output.stderr));

    // Gather blast results
    let blasts: Vec<Blast> = gather_blast_results(&outfile, &search);
    info!(
        "Found {} blast results for hmm search id# {}, gene {}",
        blasts.len().to_string(),
        search.id,
        search.gene_id
    );

    // Delete tmpfile
    fs::remove_file(&tmpfile).expect("Unable to delete temporary query file");

    // Return
    blasts
}

fn generate_tmpfile(search: &HmmSearch) -> String {
    // Get filename
    let tmpfile = format!(
        "{}/tmp/{}.hmm-{}",
        CONFIG.report.output_dir, search.gene_id, search.id
    );

    // Save to tmpfile
    let path = Path::new(&tmpfile);
    let mut fh = match File::create(&path) {
        Ok(res) => res,
        Err(e) => panic!(
            "Unable to open temporary file for writing, {}, error: {}",
            tmpfile, e
        ),
    };

    // Write to file
    let line = format!(">{}\n{}\n", search.target, search.sequence);
    fh.write_all(&line.as_bytes())
        .expect("Unable to write to temporary file");

    // Return
    tmpfile
}

fn gather_blast_results(outfile: &String, search: &HmmSearch) -> Vec<Blast> {
    // Initialize
    let mut blast_res: Vec<Blast> = Vec::new();

    // Open file
    let fh = match File::open(&outfile) {
        Ok(res) => res,
        Err(e) => panic!(
            "Unable to open blast result file at {}, error: {}",
            outfile, e
        ),
    };
    let lines = io::BufReader::new(fh).lines();

    // Go through lines
    for ln in lines {
        let line = match ln {
            Ok(l) => l,
            Err(_e) => continue,
        };

        if line.starts_with("#") {
            continue;
        }

        // Set variables
        let parts = line.trim_end().split("\t").collect::<Vec<&str>>();
        let tmp_evalue = parts[3].parse::<f32>().unwrap();
        let log_evalue = if tmp_evalue == 0.00 {
            -999.00
        } else {
            tmp_evalue.log2().floor()
        };

        // Add to blast
        blast_res.push(Blast {
            query: parts[0].to_string(),
            target: parts[1].parse::<u16>().unwrap(),
            score: parts[3].parse::<f32>().unwrap(),
            evalue: parts[2].parse::<f32>().unwrap(),
            log_evalue: log_evalue,
            res_start: parts[4].parse::<u16>().unwrap(),
            res_end: parts[5].parse::<u16>().unwrap(),
            hmmsearch_id: search.id,
        });
    }

    // Return
    blast_res
}
