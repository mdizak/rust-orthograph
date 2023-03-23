use crate::reporter::ReporterKit;
use crate::stats::Stats;
use biotools::db::sqlite::TBL_HITS;
use biotools::CONFIG;
use log::{info, warn};
use rusqlite::{Error, Statement};
use std::collections::HashMap;

struct HmmCandidate {
    id: u32,
    hmmsearch_id: u32,
    gene_id: String,
    score: f32,
    hmm_start: u16,
    hmm_end: u16,
    header_base: String,
    header_revcomp: u8,
    header_translate: u8,
    rank: u16,
}

pub struct HmmDiscard {
    pub hit_id: u32,
    pub gene_id: String,
    pub header_base: String,
    pub revcomp: u8,
    pub translate: u8,
}

pub fn check(kit: &ReporterKit, mut stats: &mut Stats) {
    // Gather candidates
    let candidates = match gather_candidates(&kit) {
        Ok(r) => r,
        Err(e) => panic!("Unable to gather hmm overlap candidates, error: {}", e),
    };

    // Process candidates
    let discards: Vec<Vec<HmmDiscard>> = candidates
        .iter()
        .map(|c| {
            let res: Vec<HmmDiscard> = process_candidates(&c.1);
            res
        })
        .collect();

    // Discard candidates
    for group in discards {
        for hit in group {
            stats.discard_hmm_overlap(&kit, &hit);
        }
    }
}

fn gather_candidates(kit: &ReporterKit) -> Result<HashMap<String, Vec<HmmCandidate>>, Error> {
    // Execute sql
    let mut stmt = prepare_sql(&kit);
    let mut rows = match stmt.query([]) {
        Ok(r) => r,
        Err(e) => panic!("Unable to execute sql for hmm overlap check, error: {}", e),
    };

    // Go through rows
    let mut candidates: HashMap<String, Vec<HmmCandidate>> = HashMap::new();
    loop {
        // Get row
        let row = match rows.next()? {
            Some(r) => r,
            None => break,
        };

        // Set hmm candidate
        let cand = HmmCandidate {
            id: row.get(0)?,
            hmmsearch_id: row.get(1)?,
            gene_id: row.get(2)?,
            score: row.get(3)?,
            hmm_start: row.get(4)?,
            hmm_end: row.get(5)?,
            header_base: row.get(6)?,
            header_revcomp: row.get(7)?,
            header_translate: row.get(8)?,
            rank: row.get(9)?,
        };

        // Add to results
        let gene_id: String = row.get(2)?;
        candidates.entry(gene_id).or_insert(Vec::new()).push(cand);
    }
    info!("Gathered all necessary candidates for hmm overlap check.");

    Ok(candidates)
}

fn prepare_sql(kit: &ReporterKit) -> Statement {
    // Set sql
    let sql = format!("
        WITH top_headers AS (
        SELECT id,hmmsearch_id,gene_id,score,hmm_start,hmm_end,header_base,header_revcomp,header_translate, RANK() OVER(PARTITION BY gene_id ORDER BY score DESC) rank FROM {}) 
        SELECT id,hmmsearch_id,gene_id,score,hmm_start,hmm_end,header_base,header_revcomp,header_translate,rank FROM top_headers ORDER BY rank
    ", *TBL_HITS);

    // Prepare sql
    let mut stmt = match kit.memdb.prepare(&sql) {
        Ok(res) => res,
        Err(e) => panic!(
            "Unable to prepare SQL statement while fetching env dupe check sequences, error: {}",
            e
        ),
    };

    stmt
}

fn process_candidates(candidates: &Vec<HmmCandidate>) -> Vec<HmmDiscard> {
    // Initialize
    let mut discards: Vec<u32> = Vec::new();
    let mut hmm_discards: Vec<HmmDiscard> = Vec::new();
    let mut hits: Vec<&HmmCandidate> = candidates.iter().rev().map(|c| c).collect();

    // Go through hits
    loop {
        // Get next hit
        let hit_a = match hits.pop() {
            Some(r) => r,
            None => break,
        };

        // Skip, if needed
        if discards.contains(&hit_a.id) {
            continue;
        }

        // Go through remaining hits
        for hit_b in &hits {
            // Skip, if needed
            if discards.contains(&hit_b.id) || hit_a.header_base == hit_b.header_base {
                continue;
            }

            // Get percent overlap
            let percent = match biotools::get_overlap_percent(
                hit_a.hmm_start..hit_a.hmm_end + 1,
                hit_b.hmm_start..hit_b.hmm_end + 1,
                false,
            ) {
                Some(r) => r,
                None => continue,
            };

            // Check overlap percent
            if percent < CONFIG.search.hmm_overlap_threshold {
                continue;
            }

            // Check score
            if (hit_a.score / hit_b.score) < CONFIG.search.hmm_score_discard_threshold {
                continue;
            }

            // Discard transcript
            hmm_discards.push(HmmDiscard {
                hit_id: hit_b.id,
                gene_id: format!("{}", hit_b.gene_id),
                header_base: format!("{}", hit_b.header_base),
                revcomp: hit_b.header_revcomp,
                translate: hit_b.header_translate,
            });
            discards.push(hit_b.id);

            warn!("Discarding hmm search {}, gene {}, header {} as it has {} percent overlap with master", hit_b.hmmsearch_id, hit_b.gene_id, hit_b.header_base, percent);
        }
    }

    hmm_discards
}
