use crate::reporter::ReporterKit;
use biotools::db::sqlite::{TBL_HITS, TBL_AASEQS, TBL_ESTS};
use crate::models::{Hit, OrfTranscript};
use crate::stats::Stats;
use crate::algorithms::{orf, orf_extended};
use rusqlite::{Error, Statement, ToSql};
use biotools::CONFIG;
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use log::warn;

pub struct OrfResult {
    pub hit_id: u32,
    orf: Option<OrfTranscript>,
    taxid: Option<u16>,
    pub gene_id: Option<String>,
    pub header_base: Option<String>,
    pub revcomp: Option<u8>,
    pub translate: Option<u8>
}

pub fn run(kit: &ReporterKit, mut stats: &mut Stats) -> Result<bool, Error> { 

    // Prepare insert sql
    let mut insert_stmt = prepare_insert_sql(&kit);

    // Prepare and execute sql
    let mut stmt = prepare_select_sql(&kit);
    let mut rows = stmt.query([])
        .expect("Unable to execute SQL to retrieve hits during finalization.");

    // GO through rows
    let mut hits: Vec<Hit> = Vec::new();
    loop {

        // Get row
        let row = match rows.next()? {
            Some(r) => r,
            None => break
        };

        // Get est sequence, reverse if necessary
        let is_revcomp: bool = row.get(21)?;
        let mut est_sequence: String = row.get(26)?;
        if is_revcomp == true {
            est_sequence = reverse_seq(&est_sequence);
        }
        let (ali_start, ali_end): (u16, u16) = (row.get(13)?, row.get(14)?);

        // Define hit
        let hit = Hit {
            id: row.get(0)?,
            is_overlap: row.get(1)?,
            hmmsearch_id: row.get(2)?,
            taxid: row.get(24)?,
            aaseq_id: row.get(4)?,
            ntseq_id: row.get(5)?,
            blast_target: row.get(6)?,
            gene_id: row.get(7)?,
            score: row.get(8)?,
            digest: row.get(9)?,
            evalue: row.get(10)?,
            hmm_start: row.get(11)?,
            hmm_end: row.get(12)?,
            ali_start: ali_start,
            ali_end: ali_end,
            env_start: row.get(15)?,
            env_end: row.get(16)?,
            blast_start: row.get(17)?,
            blast_end: row.get(18)?,
            header_base: row.get(19)?,
            header_full: row.get(20)?,
            header_revcomp: is_revcomp,
            header_translate: row.get(22)?, 
            non_orf_sequence: row.get(23)?,
            aa_sequence: row.get(25)?,
            hmm_sequence: est_to_hmm(&est_sequence, &ali_start, &ali_end),
            est_sequence: est_sequence
        };

        // Add to hits
        hits.push(hit);
    }

    // Process candidates
    let results: Vec<OrfResult> = hits.par_iter().map(|h| {
        let res = process_hits(&h);
        res
    }).collect();

    // Go through results
    for res in results {

        // Check for discard
        if None == res.taxid {
            stats.discard_non_orf(&kit, &res);
            continue;
        }
        let orf = res.orf.unwrap();

        // Insert into temp orf table
        match insert_stmt.execute([
            &orf.hit_id as &dyn ToSql,
            &res.taxid.unwrap(),
            &orf.cdna_start,
            &orf.cdna_end,
            &orf.aa_start,
            &orf.aa_end,
            &orf.cdna_start_transcript,
            &orf.cdna_end_transcript,
            &orf.aa_start_transcript,
            &orf.aa_end_transcript,
            &orf.aa_start_hmm,
            &orf.aa_end_hmm,
            &orf.translated_seq,
            &orf.cdna_seq
        ]) {
            Ok(r) => r,
            Err(e) => panic!("Unable to insert into temporary orf table, error: {}", e)
        };
    }

    Ok(true)
}

fn process_hits(hit: &Hit) -> OrfResult {

    // Generate orf
    let initial_orf = match orf::generate(&hit, false) {
        Some(r) => r,
        None => {
            warn!("Unable to generate orf for hmm search id# {}, hdr {}, gene {}, skipping transcript.", hit.hmmsearch_id, hit.header_base, hit.gene_id);
            return OrfResult {
                hit_id: hit.id,
                orf: None,
                taxid: None,
                gene_id: Some(format!("{}", hit.gene_id)),
                header_base: Some(format!("{}", hit.header_base)),
                revcomp: Some(hit.header_revcomp as u8),
                translate: Some(hit.header_translate)
            };
        }
    };

    // Check for extended orf
    let orf = match orf_extended::generate(&hit, &initial_orf) {
        Some(r) => r,
        None => {
            warn!("Did not receive extended orf, reverting to initial orf");
            initial_orf
        }
    };

    // Define orf
    let res = OrfResult {
        hit_id: hit.id,
        orf: Some(orf),
        taxid: Some(hit.taxid),
        gene_id: Some(format!("{}", hit.gene_id)),
        header_base: Some(format!("{}", hit.header_base)),
        revcomp: Some(hit.header_revcomp as u8),
        translate: Some(hit.header_translate)
    };

    res
}

fn prepare_select_sql(kit: &ReporterKit) -> Statement {

    let sql = format!("SELECT 
        h.*,
        a.taxid,
        a.sequence aa_sequence, 
        e.sequence 
    FROM 
        {} h, {} a, out.{} e 
    WHERE 
        h.blast_target = a.id AND 
        e.header = h.header_base
        ORDER BY h.id 
    ", *TBL_HITS, *TBL_AASEQS, *TBL_ESTS);
//ORDER BY h.gene_id,h.score,h.header_base
    // Prepare sql
    let stmt = match kit.memdb.prepare(&sql) {
        Ok(r) => r,
        Err(e) => panic!("Unable to prepare sql to select hits during frameshift correction, error: {}", e)
    };

    stmt
}

fn prepare_insert_sql(kit: &ReporterKit) -> Statement {

    // Format sql
    let sql = format!("INSERT INTO {}_orf (hit_id, taxid, cdna_start, cdna_end, aa_start, aa_end, cdna_start_transcript, cdna_end_transcript, aa_start_transcript, aa_end_transcript, aa_start_hmm, aa_end_hmm, translated_seq, cdna_seq) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", CONFIG.db.table_prefix);

    // Prepare
    let stmt = match kit.memdb.prepare(&sql) {
        Ok(r) => r,
        Err(e) => panic!("Unable to prepare insert sql statement for orf table, error: {}", e)
    };

    stmt
}


fn reverse_seq(old_seq: &String) -> String {

    let seq: String = old_seq.chars().rev().map(|c| match c {
        'A' => 'T',
        'G' => 'C',
        'C' => 'G',
        'T' => 'A',
        'Y' => 'R',
        'R' => 'Y',
        'K' => 'M',
        'M' => 'K',
        'a' => 'T',
        'g' => 'C',
        'c' => 'G',
        't' => 'A',
        _ => c
    }).collect();

    seq
}

fn est_to_hmm(est_sequence: &String, ali_start: &u16, ali_end: &u16) -> String {

    // Get start and end
    let start: usize = (*ali_start as usize - 1) * 3;
    let end: usize = start + ((*ali_end as usize - *ali_start as usize + 1) * 3);

    // Return
    est_sequence.as_str()[start..end].to_string()
}


