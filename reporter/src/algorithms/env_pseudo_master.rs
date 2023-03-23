use crate::reporter::ReporterKit;
use crate::stats::Stats;
use biotools::db::sqlite::TBL_HITS;
use rusqlite::{Statement, Error};
use std::collections::HashMap;
use biotools::CONFIG;
use log::{info, warn};

pub struct EnvCandidate {
    pub id : u32,
    hmmsearch_id: u32,
    pub gene_id: String,
    score: f32,
    env_start: u16,
    env_end: u16,
    pub header_base: String,
    pub hdr_revcomp: u8,
    pub hdr_translate: u8,
    rank: u8
}

pub fn check(kit: &ReporterKit, mut stats: &mut Stats) -> Result<bool, Error> {

    // Execute sql
        let mut stmt = prepare_sql(&kit);
    let mut rows = match stmt.query([]) {
        Ok(res) => res,
        Err(e) => panic!("Unable to execute sql while fetching filtered scores from database, error: {}", e)
    };
    let mut hits: HashMap<String, Vec<EnvCandidate>> = HashMap::new();

    // Go through rows
    loop {

        // Get next row
        let row = match rows.next()? {
            Some(r) => r,
            None => break
        };

        // Define hit
        let hit = EnvCandidate {
            id: row.get(0)?,
            hmmsearch_id: row.get(1)?,
            gene_id: row.get(2)?,
            score: row.get(3)?,
            env_start: row.get(4)?,
            env_end: row.get(5)?,
            header_base: row.get(6)?,
            hdr_revcomp: row.get(7)?,
            hdr_translate: row.get(8)?,
            rank: row.get(9)?
        };

        // Add to hits
        let header = format!("{}", hit.header_base);
        hits.entry(header).or_insert(Vec::new()).push(hit);
    }

    // Process
    let results: Vec<bool> = hits.iter().map(|c| {
        let res = process_candidates(&kit, &c.1, &mut stats);
        res
    }).collect();

    Ok(true)
}

fn process_candidates(kit: &ReporterKit, candidates: &Vec<EnvCandidate>, mut stats: &mut Stats) -> bool {

    // Get master
    let master = candidates.first().unwrap();
    let mut master_start: u16 = master.env_start;
    let mut master_end: u16 = master.env_end;
    let mut is_minescule: bool = false;
    info!("Checking master-pseudo for base header {}, gene {} which has {} child transcripts with the same base header.", master.header_base, master.gene_id, candidates.len() - 1);

    // Extend master coords, if necessary
    match extend_master_coords(&kit, &master, &candidates) {
        Some(r) => {
            master_start = r.0;
            master_end = r.1;
        },
        None => { }
    };

    // GO through non-gene sequences
    for cand in candidates {

        // Skip, if master
        if cand.id == master.id {
            continue;
        }

        // Discard, if same gene as master
        if cand.gene_id == master.gene_id {
            stats.discard_env_pseudo_master(&kit, &cand);
            continue;
        }

        // Get overlap percent
        let percent = match biotools::get_overlap_percent(master_start..master_end + 1, cand.env_start..cand.env_end + 1, true) {
            Some(r) => r,
            None => continue
        };

        // Check percent
        if percent < CONFIG.search.env_overlap_threshold {
            is_minescule = true;
            break;
        }

            // Check score
            if master.score / cand.score >= CONFIG.search.env_score_discard_threshold {
            warn!("child transcript of base header {} in gene {} has overlap of {} and score of {}, discarding transcript.", cand.header_base, cand.gene_id, percent, cand.score);
            stats.discard_env_overlap(&kit, &cand);
        } else {
            info!("Transcript hdr {} in gene {} only overlaps master by {} percent, keeping transcript.", cand.header_base, cand.gene_id, percent);
        }

    }

    // Return if needed
    if !is_minescule {
        return true;
    }

    // Go through candidates
    for cand in candidates {
        stats.discard_env_overlap(&kit, &cand);
    }

    true
}

fn extend_master_coords(kit: &ReporterKit, master: &EnvCandidate, children: &Vec<EnvCandidate>) -> Option<(u16, u16)> {

    // Get pseudo masters
    let mut mst_start: u16 = master.env_start;
    let mut mst_end: u16 = master.env_end;

    // Go through children
    for c in children {

        if c.env_start < mst_start {
            mst_start = c.env_start;
        }

        if c.env_end > mst_start {
            mst_end = c.env_end;
        }
    }

    // Check for no change
    if mst_start == master.env_start && mst_end == master.env_end {
        return None;
    }

    // Update master env_start coord
    if mst_start > 0 && mst_start < master.env_start {
        let sql = format!("UPDATE {} SET env_start = {} WHERE id = ?", *TBL_HITS, mst_start);
        if let Err(e) = kit.memdb.execute(&sql, [&master.id]) {
            panic!("Unable to update env_start on master during env overlap check, error: {}", e);
        }
    }

    // Update master env_start coord
    if mst_end > master.env_end {
        let sql = format!("UPDATE {} SET env_end = {} WHERE id = ?", *TBL_HITS, mst_end);
        if let Err(e) = kit.memdb.execute(&sql, [&master.id]) {
            panic!("Unable to update env_end on master during env overlap check, error: {}", e);
        }
    }

    Some((mst_start, mst_end))
}

fn prepare_sql(kit: &ReporterKit) -> Statement {

    // Set sql
    let sql = format!("
        WITH top_headers AS (
        SELECT id,hmmsearch_id,gene_id,score,env_start,env_end,header_base,header_revcomp,header_translate, RANK() OVER(PARTITION BY header_base ORDER BY score DESC) rank FROM {}) 
        SELECT id,hmmsearch_id,gene_id,score,env_start,env_end,header_base,header_revcomp,header_translate,rank FROM top_headers ORDER BY header_base,rank
    ", *TBL_HITS);

    // Prepare sql
    let stmt = match kit.memdb.prepare(&sql) {
        Ok(res) => res,
        Err(e) => panic!("Unable to prepare SQL statement while fetching env dupe check sequences, error: {}", e)
    };

    stmt
}

