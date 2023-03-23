use crate::models::{Hit, OrfTranscript};
use biotools::{FastaResult, CONFIG};
use log::{info, warn};
use std::collections::HashMap;
use std::process::Command;
use std::string::String;

pub fn generate(hit: &Hit, is_extended: bool) -> Option<OrfTranscript> {
    // Check for no framework-correction
    if !CONFIG.switch.frameshift_correction {
        return Some(OrfTranscript {
            hit_id: hit.id,
            translated_seq: format!("{}", hit.non_orf_sequence),
            cdna_seq: format!("{}", hit.hmm_sequence),
            cdna_start: (hit.ali_start - 1) * 3,
            cdna_end: ((hit.ali_start - 1) * 3) + ((hit.ali_end - hit.ali_start + 1) * 3),
            aa_start: hit.ali_start,
            aa_end: hit.ali_end,
            cdna_start_transcript: (hit.ali_start - 1) * 3,
            cdna_end_transcript: ((hit.ali_start - 1) * 3)
                + ((hit.ali_end - hit.ali_start + 1) * 3),
            aa_start_transcript: hit.ali_start,
            aa_end_transcript: hit.ali_end,
            aa_start_hmm: hit.hmm_start,
            aa_end_hmm: hit.hmm_end,
        });
    }

    // Save query file
    let query_contents = format!(">query\n{}", hit.aa_sequence);
    let query_file = biotools::io::create_tmp_file(&query_contents);

    // Save target file
    let target_contents = if is_extended {
        format!(">target\n{}\n", hit.est_sequence)
    } else {
        format!(">target\n{}\n", hit.hmm_sequence)
    };
    let target_file = biotools::io::create_tmp_file(&target_contents);

    // Run exxonerate command
    let output = Command::new(&CONFIG.programs.exonerate)
        .args([
            "--bestn",
            "1",
            "--score",
            &CONFIG.search.hmmsearch_threshold.to_string(),
            "--ryo",
            ">cdna %tcb %tce\n%tcs>aa %qab %qae\n%qas",
            "--subopt",
            "0",
            "--geneticcode",
            "1",
            "--model",
            "protein2genome",
            "--querytype",
            "protein",
            "--targettype",
            "dna",
            "--verbose",
            "0",
            "--showalignment",
            "no",
            "--showvulgar",
            "no",
            "--query",
            &query_file,
            "--target",
            &target_file,
        ])
        .output()
        .expect("Unable to execute exonerate program.");

    // Check output status
    if !output.status.success() {
        warn!(
            "Did not receive successful exit code from exonerate, skipping transcript.  Error: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        return None;
    }

    // Log
    info!(
        "Completed exonerate search of hmm search id# {} gene {} with {}, (stderr: {})",
        hit.hmmsearch_id,
        hit.gene_id,
        output.status,
        String::from_utf8_lossy(&output.stderr)
    );

    // Read output
    let output_contents = String::from_utf8_lossy(&output.stdout);
    let fasta = biotools::read_fasta_string(&output_contents.to_string());

    // Read exonerate response
    let (cdna, aa) = match read_exonerate_response(&fasta) {
        Some(r) => r,
        None => return None,
    };

    // Translate
    let translated = match translate(&cdna.sequence) {
        Some(r) => r,
        None => {
            warn!("Did not receive a valid response from translate program, skipping transcript.");
            return None;
        }
    };

    // Return
    if is_extended {
        return Some(OrfTranscript {
            hit_id: hit.id,
            translated_seq: translated,
            cdna_seq: format!("{}", cdna.sequence),
            cdna_start: cdna.coord_start + 1,
            cdna_end: cdna.coord_end,
            aa_start: aa.coord_start,
            aa_end: aa.coord_end,
            cdna_start_transcript: cdna.coord_start + 1 + (hit.ali_start * 3) - 3,
            cdna_end_transcript: cdna.coord_end + (hit.ali_start * 3) - 3,
            aa_start_transcript: (cdna.coord_start / 3 + 1),
            aa_end_transcript: (cdna.coord_end / 3) + 1,
            aa_start_hmm: hit.hmm_start + (cdna.coord_start + ((hit.ali_start * 3) - 3) / 3),
            aa_end_hmm: hit.hmm_start + ((cdna.coord_end + ((hit.ali_start * 3) - 3)) / 3),
        });
    }

    Some(OrfTranscript {
        hit_id: hit.id,
        translated_seq: translated,
        cdna_seq: format!("{}", cdna.sequence),
        cdna_start: cdna.coord_start + 1,
        cdna_end: cdna.coord_end,
        aa_start: aa.coord_start,
        aa_end: aa.coord_end,
        cdna_start_transcript: cdna.coord_start + 1 + (hit.ali_start * 3) - 3,
        cdna_end_transcript: cdna.coord_end + (hit.ali_start * 3) - 3,
        aa_start_transcript: ((cdna.coord_start + ((hit.ali_start * 3) - 3) / 3) as f32).floor()
            as u16,
        aa_end_transcript: (((cdna.coord_end + (((hit.ali_start + 1) * 3) - 3)) / 3) as f32).ceil()
            as u16,
        aa_start_hmm: hit.hmm_start + (cdna.coord_start + ((hit.ali_start * 3) - 3) / 3),
        aa_end_hmm: hit.hmm_start + ((cdna.coord_end + ((hit.ali_start * 3) - 3)) / 3),
    })
}

fn read_exonerate_response<'a>(
    fasta: &'a HashMap<String, FastaResult>,
) -> Option<(&'a FastaResult, &'a FastaResult)> {
    // Get initial cdna
    let cdna = match fasta.get(&"cdna1".to_string()) {
        Some(r) => r,
        None => return None,
    };

    // Get initial aa
    let aa = match fasta.get(&"aa2".to_string()) {
        Some(r) => r,
        None => return None,
    };

    Some((cdna, aa))
}

fn translate(cdna_seq: &String) -> Option<String> {
    // Save tmp file
    let tmp_contents = format!(">cdna\n{}", cdna_seq);
    let tmpfile = biotools::io::create_tmp_file(&tmp_contents);

    // Execute fasta-translate command
    let output = Command::new(&CONFIG.programs.translate)
        .args(["--geneticcode", "1", "-F", "1", &tmpfile])
        .output()
        .expect("Unable to run translate program");

    // Check status
    if !output.status.success() {
        warn!("Did not receive successful exit status from translate program.  Got {} with (stdout: {}), (stderr: {})", 
            output.status, String::from_utf8_lossy(&output.stdout), String::from_utf8_lossy(&output.stderr));
        return None;
    }

    // Read fasta file
    let output_contents = String::from_utf8_lossy(&output.stdout);
    let fasta = biotools::read_fasta_string(&output_contents.to_string());
    let cdna = match fasta.get(&"cdna1".to_string()) {
        Some(r) => r,
        None => return None,
    };

    let res_seq = format!("{}", &cdna.sequence);
    Some(res_seq)
}
