use rusqlite::{Connection, Error};
use std::path::Path;
use std::collections::HashMap;
use lazy_static::lazy_static;
use crate::db::models::BlastResult;
use crate::CONFIG;
use log::debug;

lazy_static! {
    pub static ref TBL_SPECIES_INFO: String = format!("{}_species_info", CONFIG.db.table_prefix);
    pub static ref TBL_ESTS: String = format!("{}_ests", CONFIG.db.table_prefix);
    pub static ref TBL_HMMSEARCH: String = format!("{}_hmmsearch", CONFIG.db.table_prefix);
    pub static ref TBL_BLAST: String = format!("{}_blast", CONFIG.db.table_prefix);
    pub static ref TBL_AASEQS: String = format!("input.{}_aaseqs", CONFIG.db.table_prefix);
    pub static ref TBL_SEQUENCE_PAIRS: String = format!("input.{}_sequence_pairs", CONFIG.db.table_prefix);
    pub static ref TBL_BLASTDBS: String = format!("input.{}_blastdbs", CONFIG.db.table_prefix);
    pub static ref TBL_SEQUENCE_TYPES: String = format!("input.{}_sequence_types", CONFIG.db.table_prefix);
    pub static ref TBL_NTSEQS: String = format!("input.{}_ntseqs", CONFIG.db.table_prefix);
    pub static ref TBL_SET_DETAILS: String = format!("input.{}_set_details", CONFIG.db.table_prefix);
    pub static ref TBL_OGS: String = format!("input.{}_ogs", CONFIG.db.table_prefix);
    pub static ref TBL_TAXA: String = format!("input.{}_taxa", CONFIG.db.table_prefix);
    pub static ref TBL_LOGS: String = format!("input.{}_orthologs", CONFIG.db.table_prefix);
    pub static ref TBL_HITS: String = format!("{}_hits", CONFIG.db.table_prefix);
}

pub struct Sqlite {
    pub conn: Connection
}

struct RefTaxa {
    taxa: String
}

struct AaseqByGene {
    gene: String,
    seqid: u32
}

impl Sqlite {

    pub fn new() -> Self {

        // Check sqlite file exists
        if !Path::new(&CONFIG.db.reporter_sqlite_file).exists() {
            panic!("SQLite database file does not exist at {}.", CONFIG.db.reporter_sqlite_file);
        }

        // Connect to database
        let conn = Connection::open(&CONFIG.db.reporter_sqlite_file).unwrap_or_else(|_error| {
            panic!("Unable to open SQLite database at {}", &CONFIG.db.reporter_sqlite_file)
        });

        // Attach input database
        let sql = format!("ATTACH '{}' AS input", CONFIG.db.sqlite_file);
        conn.execute(&sql, [])
            .expect("Unable to attach input SQLite database");

        // Return
        Self {
            conn: conn
        }

    }

    pub fn get_species_id(&self) -> Result<u32, rusqlite::Error> { 

        // Prepare
        let sql = format!("SELECT id,name FROM {} WHERE name = ?", *TBL_SPECIES_INFO);
        let mut stmt = match self.conn.prepare(&sql) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        // Execute
        let mut rows = match stmt.query([&CONFIG.report.species_name]) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        // Get row
        let row = match rows.next() {
            Ok(res) => res.unwrap(),
            Err(error) => return Err(error)
        };

        // Get id#
        let id: u32 = match row.get(0) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        Ok(id)
    }

    pub fn get_set_id(&self) -> Result<u32, rusqlite::Error> { 

        // Prepare
        let sql = format!("SELECT id,name FROM {} WHERE name = ?", *TBL_SET_DETAILS);
        let mut stmt = match self.conn.prepare(&sql) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        // Execute
        let mut rows = match stmt.query([&CONFIG.report.set_name]) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        // Get row
        let row = match rows.next() {
            Ok(res) => res.unwrap(),
            Err(error) => return Err(error)
        };

        // Get id#
        let id: u32 = match row.get(0) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        Ok(id)
    }

    pub fn get_reference_taxa(&self, set_id: &u32) -> Result<Vec<String>, Error> {

        // Set SQL
        let sql = format!("SELECT DISTINCT t.name
            FROM {} p, {} t, {} l 
            WHERE t.id = p.taxid AND l.sequence_pair = p.id AND l.setid = ?", 
        *TBL_SEQUENCE_PAIRS, *TBL_TAXA, *TBL_LOGS);

        // Prepare
        let mut stmt = match self.conn.prepare(&sql) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        // Execute
        let rows = stmt.query_map([&set_id], |row| {
            Ok( RefTaxa { taxa: row.get(0)? })
        })?;

        // Collect, and return
        let taxa: Vec<String> = rows.map(|r|r.unwrap().taxa).collect::<Vec<String>>();
        Ok(taxa)
    }

    pub fn get_aaseq_in_set(&self, set_id: &u32) -> Result<HashMap<String, Vec<u32>>, Error> {

        // Format sql
        let sql = format!("SELECT DISTINCT l.ortholog_gene_id, a.id  
            FROM {} l, {} a, {} p 
            WHERE l.sequence_pair = p.id AND p.aa_seq = a.id AND l.setid = ?", 
        *TBL_LOGS, *TBL_AASEQS, *TBL_SEQUENCE_PAIRS);

        // Prepare
        let mut stmt = match self.conn.prepare(&sql) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        // Execute
        let rows = stmt.query_map([&set_id], |row| {
            Ok( AaseqByGene { 
                gene: row.get(0)?,
                seqid: row.get(1)?
            })
        })?;

        // Collect, and return
        let mut aaseq: HashMap<String, Vec<u32>> = HashMap::new();
        for r in rows {
            let row = r.unwrap();
            aaseq.entry(row.gene).or_insert(Vec::new()).push(row.seqid);
        }

        // Return
        Ok(aaseq)
    }

    pub fn get_hmm_blast_results(&self, hmmsearch_id: u32) -> Result<Vec<BlastResult>, Error> {

        // Set sql
        let sql = format!("SELECT DISTINCT 
            b.target,
            b.score,
            b.evalue,
            b.start,
            b.end
            FROM {} b, {} s, {} e 
            WHERE  
            s.id = b.hmmsearch_id AND 
            e.digest = s.target AND 
            s.target IS NOT NULL AND 
            b.hmmsearch_id = ? 
            ORDER BY b.score DESC",
        *TBL_BLAST, *TBL_HMMSEARCH, *TBL_ESTS);

        // Prepare
        let mut stmt = match self.conn.prepare(&sql) {
            Ok(res) => res,
            Err(error) => return Err(error)
        };

        // Execute
        let rows = stmt.query_map([&hmmsearch_id], |row| {
            Ok( BlastResult { 
                target: row.get(0)?,
                score: row.get(1)?,
                evalue: row.get(2)?,
                res_start: row.get(3)?,
                res_end: row.get(4)?
            })
        })?;

        let res: Vec<BlastResult> = rows.into_iter().map(|r| r.unwrap()).collect();

        // Return
        Ok(res)
    }

    pub fn get_ref_taxon_name(&self, aaseq_id: &u32) -> Result<String, Error> {

        // Set sql
        let sql = format!("SELECT t.name FROM {} t, {} a WHERE t.id = a.taxid AND a.id = ?", *TBL_TAXA, *TBL_AASEQS);

        // Prepare
        let mut stmt = match self.conn.prepare(&sql) {
            Ok(res) => res,
            Err(e) => panic!("Unable to prepare sql statement to retrieve reference taxon name, error: {}", e)
        };

        // Execute
        let mut rows = match stmt.query([&aaseq_id]) {
            Ok(res) => res,
            Err(e) => panic!("Unable to retrive ref taxon name of aaseq id# {}, error: {}", &aaseq_id, e)
        };

        // Get row
        let row = match rows.next() {
            Ok(res) => res.unwrap(),
            Err(e) => panic!("Unable to retrive ref taxon name of aaseq id# {}, error: {}", &aaseq_id, e)
        };

        // Return
        let name = row.get(0)?;
        debug!("Obtained ref taxon name {} for aaseq id# {}", &name, &aaseq_id);
        Ok(name)
    }

    pub fn get_blast_count(&self, search_id: &u32) -> u16 {

        // Prepare sql
        let sql = format!("SELECT count(*) FROM {} WHERE hmmsearch_id = ?", *TBL_BLAST);
        let mut stmt = match self.conn.prepare(&sql) {
            Ok(res) => res,
            Err(e) => panic!("Unable to prepare sql to get blast result count, error: {}", e)
        };

        // Execute
        let mut rows = match stmt.query([&search_id]) {
            Ok(r) => r,
            Err(e) => panic!("Unable to execute sql to get blast count, error: {}", e)
        };

        // Get row
        let row = rows.next().unwrap().unwrap();
        row.get(0).unwrap()
    }


}


