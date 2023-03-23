use crate::models::Blast;
use crate::BIODB_ARGS;
use rusqlite::{Connection, Error};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

pub struct Sqlite {
    pub conn: Arc<Connection>,
}

impl Sqlite {
    pub fn new() -> Self {
        // Check sqlite file exists
        if !Path::new(&BIODB_ARGS.sqlite_file).exists() {
            panic!(
                "SQLite database file does not exist at {}.",
                BIODB_ARGS.sqlite_file
            );
        }

        // Connect to database
        let conn = Connection::open(&BIODB_ARGS.sqlite_file).unwrap_or_else(|_error| {
            panic!(
                "Unable to open SQLite database at {}",
                &BIODB_ARGS.sqlite_file
            )
        });

        // Return
        Self {
            conn: Arc::new(conn),
        }
    }

    pub fn get_total(&self, table_name: &str) -> u64 {
        let sql = format!("SELECT count(*) FROM {}", table_name);
        let total: u64 = self
            .conn
            .query_row::<u64, _, _>(&sql, [], |row| row.get(0))
            .unwrap();
        total
    }

    pub fn get_blast_results(&self, hmmsearch_id: &u32) -> Result<Vec<Blast>, Error> {
        // Execute sql
        let mut stmt = self.conn.prepare("SELECT b.*,e.header,e.type FROM orthograph_blast b, orthograph_ests e WHERE b.hmmsearch_id = ? AND b.query = e.digest ORDER BY score DESC").unwrap();
        let mut rows = match stmt.query([&hmmsearch_id]) {
            Ok(r) => r,
            Err(e) => panic!(
                "Unable to execute SQL to obtain blast results, error: {}",
                e
            ),
        };
        let mut blast_results: Vec<Blast> = Vec::new();

        // Go through rows
        loop {
            // Get row
            let row = match rows.next()? {
                Some(r) => r,
                None => break,
            };
            let header: String = row.get(10)?;

            let blast = Blast {
                taxid: row.get(1)?,
                target: row.get(3)?,
                score: row.get(4)?,
                evalue: row.get(5)?,
                log_evalue: row.get(6)?,
                blast_start: row.get(7)?,
                blast_end: row.get(8)?,
                seq_type: row.get(11)?,
                header: header.replace(" ", "_").to_string(),
            };
            blast_results.push(blast);
        }

        Ok(blast_results)
    }
}
