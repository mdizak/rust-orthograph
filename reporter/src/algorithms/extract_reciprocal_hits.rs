use crate::reporter::ReporterKit;
use biotools::db::sqlite::{TBL_HMMSEARCH, TBL_LOGS, TBL_ESTS, TBL_HITS, TBL_SEQUENCE_PAIRS};
use crate::models::HmmSearch;
use rusqlite::{Statement, Error};
use crate::algorithms::is_reciprocal_hit;
use crate::stats::Stats;
use biotools::CONFIG;
use log::{info, warn};
use rusqlite::ToSql;


pub fn run(kit: &ReporterKit) -> Result<Stats, Error> {

    // Prepare and execute sql
    let mut stmt = prepare_sql(&kit);
    let mut rows = match stmt.query([&kit.species_id, &kit.set_id]) {
        Ok(res) => res,
        Err(e) => panic!("Unable to execute sql while fetching filtered scores from database, error: {}", e)
    };

    // Prepare insert sql statement
    let insert_sql = format!("INSERT INTO {} (hmmsearch_id, taxid, aaseq_id, ntseq_id, blast_target, gene_id, score, digest, evalue, hmm_start, hmm_end, ali_start, ali_end, env_start, env_end, blast_start, blast_end, header_base, header_full, header_revcomp, header_translate, non_orf_sequence) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", *TBL_HITS);
    let mut insert_stmt = match kit.memdb.prepare(&insert_sql) {
        Ok(res) => res,
        Err(e) => panic!("Unable to prepare SQL statement for insert into hits table, error: {}", e)
    };

    // Instantiate stats
    let mut stats = Stats::new();

    // Go through rows
    loop {

        // Get next row
        let row = match rows.next()? {
            Some(r) => r,
            None => break
        };

        // Define search
        let cand = HmmSearch {
            gene_id: row.get(0)?,
            aaseq_id: row.get(1)?,
            ntseq_id: row.get(2)?,
            taxid: row.get(3)?,
            hmm_id: row.get(4)?,
            score: row.get(5)?,
            digest: row.get(6)?,
            evalue: row.get(7)?,
            hmm_start: row.get(8)?,
            hmm_end: row.get(9)?,
            ali_start: row.get(11)?,
            ali_end: row.get(12)?,
            env_start: row.get(13)?,
            env_end: row.get(14)?,
            header: row.get(15)?,
            non_orf_sequence: row.get(16)?
        };

        // Skip, if not in list of wanted genes
        if CONFIG.report.wanted_genes.len() > 0 && !CONFIG.report.wanted_genes.contains(&cand.gene_id) {
            warn!("Not in list of wanted genes, skipping {}", cand.gene_id);
            continue;
        }

        // Check if reciprocal hit
        let blast = match is_reciprocal_hit::check(&kit, &cand) {
            Some(r) => r,
            None => {
                warn!("No orthology detected for {}.", cand.gene_id);
                stats.add_non_reciprocal_hit(&cand);
                continue;
            }
        };

        // Success message
        info!("Orthology detected for {}! Queueing for further checks: {}[{}:{}] to {}.", 
            cand.gene_id, cand.header, cand.hmm_start, cand.hmm_end, cand.gene_id);

        // Translate header
        let (hdr_base, hdr_revcomp, hdr_translate) = biotools::translate_header(&cand.header);

        // Insert hit to temp database table
        match insert_stmt.execute([
            &cand.hmm_id as &dyn ToSql,
            &cand.taxid,
            &cand.aaseq_id,
            &cand.ntseq_id,
            &blast.0,
            &cand.gene_id,
            &cand.score,
            &cand.digest,
            &cand.evalue,
            &cand.hmm_start,
            &cand.hmm_end,
            &cand.ali_start,
            &cand.ali_end,
            &cand.env_start,
            &cand.ali_end,
            &blast.1,
            &blast.2,
            &hdr_base.trim_end(),
            &cand.header.trim_end_matches(" "),
            &hdr_revcomp,
            &hdr_translate,
            &cand.non_orf_sequence
        ]) {
            Ok(res) => res,
            Err(e) => panic!("Unable to insert into temporary hits table with error: {}", e)
        };

    }

    // Return
    Ok(stats)
}

fn prepare_sql(kit: &ReporterKit) -> Statement {

    let sql = format!("SELECT DISTINCT 
        l.ortholog_gene_id,
        p.aa_seq,
        p.nt_seq,
        s.taxid,
        s.id,
        s.score,
        s.target,
        s.evalue,
        s.hmm_start,
        s.hmm_end,
        (s.ali_end - s.ali_start) ali_length,
        s.ali_start,
        s.ali_end,
        s.env_start,
        s.env_end,
        e.header,
        substr(e.sequence, s.ali_start, (s.ali_end - s.ali_start + 1))
        FROM {} s, {} l, {} e, {} p  
        WHERE 
            s.target = e.digest AND 
            s.query = l.ortholog_gene_id AND  
            l.sequence_pair = p.id AND 
            e.digest IS NOT NULL AND 
            s.score >= {} AND 
            s.taxid = ? AND 
            l.setid = ? AND 
            (s.ali_end - s.ali_start) + 1 >= {} 
        GROUP BY s.id ORDER BY s.score DESC", 
        *TBL_HMMSEARCH, *TBL_LOGS, *TBL_ESTS, *TBL_SEQUENCE_PAIRS, &CONFIG.search.hmmsearch_threshold, &CONFIG.search.min_transcript_length);

        // Prepare
    let stmt = match kit.db.conn.prepare(&sql) {
        Ok(res) => res,
        Err(e) => panic!("Unable to prepare SQL statement to extract reciprocal hits, error: {}", e)
    };

    // Return
    stmt
}


