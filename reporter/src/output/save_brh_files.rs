use crate::reporter::ReporterKit;
use crate::stats::Stats;
use crate::algorithms::region_mapped_before;
use crate::models::Hit;
use biotools::db::sqlite::TBL_HITS;
use std::collections::HashMap;
use std::ops::Range;
use rusqlite::Error;

pub fn save(kit: &ReporterKit, mut stats: &mut Stats) -> Result<bool, Error> {

    // Initialize
    let mut coords: HashMap<String, Vec<Range<u16>>> = HashMap::new();

    // Prepare
    let sql = format!("SELECT * FROM {} ORDER BY score DESC", *TBL_HITS);
    let mut stmt = kit.memdb.prepare(&sql)
        .expect("Unable to prepare SQL statement to retrieve all hits from temporary table.");

    // Execute sql
    let mut rows = stmt.query([])
        .expect("Unable to execute SQL to retrieve all hits from temporary table.");

    // Go through rows
    loop {

        let row = match rows.next()? {
            Some(r) => r,
            None => break
        };

        // Create hit
        let hit = Hit {
            id: row.get(0)?,
            is_overlap: row.get(1)?,
            hmmsearch_id: row.get(2)?,
            taxid: row.get(3)?,
            aaseq_id: row.get(4)?,
            ntseq_id: row.get(5)?,
            blast_target: row.get(6)?,
            gene_id: row.get(7)?,
            score: row.get(8)?,
            digest: row.get(9)?,
            evalue: row.get(10)?,
            hmm_start: row.get(11)?,
            hmm_end: row.get(12)?,
            ali_start: row.get(13)?,
            ali_end: row.get(14)?,
            env_start: row.get(15)?,
            env_end: row.get(16)?,
            blast_start: row.get(17)?,
            blast_end: row.get(18)?,
            header_base: row.get(19)?,
            header_full: row.get(20)?,
            header_revcomp: row.get(21)?,
            header_translate: row.get(22)?,
            non_orf_sequence: row.get(23)?, 
            est_sequence: "".to_string(),
            hmm_sequence: "".to_string(),
            aa_sequence: "".to_string()
        };

        // Write to brh file
        stats.write_brh(&hit);

        // Get non-overlap coords
        let digest = format!("{}", hit.digest);
        let chk_coords = coords.entry(digest).or_insert(Vec::new());

        // Add to non-overlapping file, if needed
        if !region_mapped_before::check(&kit, &hit, &chk_coords) {
            stats.write_nolap(&hit);
        }
        coords.entry(hit.digest).or_insert(Vec::new()).push(hit.ali_start..hit.ali_end);
    }

        // Return
    Ok(true)

}



