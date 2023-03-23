use biotools::db::sqlite::Sqlite;
use crate::blast;
use rusqlite::{Error, ToSql};
use crate::models::{HmmSearch, Blast};
use biotools::db::sqlite::{TBL_HMMSEARCH, TBL_ESTS, TBL_BLAST};
use rayon::iter::IntoParallelRefIterator;
use rayon::iter::ParallelIterator;
use log::{info, warn};
use crate::CONFIG;

pub struct Blastpal { 
    db: Sqlite,
    species_id: u32
}

impl Blastpal {

    pub fn new() -> Self {

        // Connect to database
    let db = Sqlite::new();
    info!("Successfully connected to SQLite database.");

        // Get species id
        let species_id: u32 = match db.get_species_id() {
            Ok(res) => res,
            Err(e) => panic!("Unable to determine id# for species {}, error: {}", &CONFIG.report.species_name, e)
        };
        info!("Got species id# {} for species name {}", species_id, CONFIG.report.species_name);

        return Self { 
            db: db,
            species_id: species_id
        }

    }

    pub fn process(self) -> Result<bool, Error> {

        // Prepare delete sql statement
        let delete_sql = format!("DELETE FROM {} WHERE hmmsearch_id = ?", *TBL_BLAST);
        let mut stmt_delete = self.db.conn.prepare(&delete_sql)
            .expect("Unable to prepare SQL to delete existing blast results");

        // Define sql 
        let hmmsearch_sql = format!("SELECT s.id,s.query,s.score,s.target,s.ali_start,s.ali_end,substr(e.sequence, s.ali_start, (s.ali_end - s.ali_start)) sequence FROM 
            {} s LEFT JOIN {} e ON e.digest = s.target AND s.score >= {} 
            LEFT JOIN {} ON s.id = hmmsearch_id GROUP BY s.id HAVING count(hmmsearch_id) NOT BETWEEN 1 AND 99", 
        *TBL_HMMSEARCH, *TBL_ESTS, CONFIG.search.hmmsearch_threshold, *TBL_BLAST);

        // Prepare SQL statement
        let mut stmt = match self.db.conn.prepare(&hmmsearch_sql) {
            Ok(r) => r,
            Err(e) => panic!("Unable to prepare sql to retrive hmmsearches, error: {}", e)
        };

        // Execute sql
        let mut rows = match stmt.query([]) {
            Ok(r) => r,
            Err(e) => panic!("Unable to execute SQL to retrieve hmmsearches, error: {}", e)
        };
        let mut pending_blasts: Vec<HmmSearch> = Vec::new();

        // GO through rows
        loop {

            // Get next row
            let row = match rows.next()? {
                Some(r) => r,
                None => break
            };

            // Create hmmsearch
            let search = HmmSearch {
                id: row.get(0)?,
                gene_id: row.get(1)?,
                target: row.get(3)?,
                ali_start: row.get(4)?,
                ali_end: row.get(5)?,
                sequence: row.get(6)?
            };
            info!("Sequence: {}", search.sequence);

            // Get total blast results
        //let total = self.db.get_blast_count(&search.id);

            // Check
            //if total != 0 && total < 100 {
                //warn!("Shkipping hmm search id# {}, gene {} as it has {} blast results, not the desired 100", search.id, search.gene_id, total);
                //continue;
            if CONFIG.report.wanted_genes.len() > 0 && !CONFIG.report.wanted_genes.contains(&search.gene_id) {
                warn!("Skipping hmm search id# {} as the gene {} is not in the wanted cog list.", search.id, search.gene_id);
                continue;
            }

            // Delete previous blast results
        stmt_delete.execute([&search.id])
                .expect("Unable to delete blast results");

            // Add to queue
            info!("Queueing hmm search id# {} for blast", search.id);
            pending_blasts.push(search);

            // Process pending blasts
            if self.run_blasts(&pending_blasts, false) {
                pending_blasts.clear();
            }

        }

        // Catch any remaining pending blasts
        self.run_blasts(&pending_blasts, true);

        // Return
        Ok(true)
    }

    fn run_blasts(&self, pending_blasts: &Vec<HmmSearch>, is_last: bool) -> bool {

        // Check length
        if pending_blasts.len() < CONFIG.search.num_threads as usize {
            return false;
        } else if pending_blasts.len() == 0 || !is_last {
            return false;
        }
        info!("Starting paralell processing of {} blasts", pending_blasts.len().to_string());

        // Process pending blasts
        let res: Vec<Vec<Blast>> = pending_blasts.par_iter().map(|s| {
            let blasts = blast::run(&s);
            blasts
        }).collect();

        // Save results
        let mut total: u16 = 0;
        for batch in res {
            total += self.save_blasts(&batch);
        } 
        info!("Completed paralell processing run, saved {} blast results to database", total);

        true
    }

    fn save_blasts(&self, blasts: &Vec<Blast>) -> u16 {

        // Prepare sql
        let insert_sql = format!("INSERT INTO {} (taxid, query, target, score, evalue, log_evalue, start, end, hmmsearch_id) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", *TBL_BLAST);
        let mut stmt = match self.db.conn.prepare(&insert_sql) {
            Ok(res) => res, 
            Err(e) => panic!("Unable to prepare SQL statement to insert into blast table, error: {}", e)
        };
        let mut total = 0;

        // GO through blasts
        for b in blasts {

            match stmt.execute([
                &self.species_id as &dyn ToSql, 
                &b.query,
                &b.target,
                &b.score,
                &b.evalue,
                &b.log_evalue,
                &b.res_start,
                &b.res_end,
                &b.hmmsearch_id
            ]) {
                Ok(res) => res,
                Err(e) => panic!("Unable to execute SQL statement to insert into blast table, error: {}", e)
            };
            total += 1;
        }

        total
    }

}

