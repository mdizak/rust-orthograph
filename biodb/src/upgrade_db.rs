extern crate serde;

use crate::models::{HmmSearch, Sequence};
use crate::sqlite::Sqlite;
use crate::database::Database;
use rusqlite::Error;
use std::io;
use std::io::Write;
use crate::ROCKSDB;
use log::info;

pub fn upgrade() {

    // Connect to SQLite
    let sqlite = Sqlite::new();

    // Transfer est sequences
    match transfer_est_sequences(&sqlite) {
        Ok(_) => { },
        Err(e) => panic!("Unable to transfer est sequences, error: {}", e)
    };

    // Transfer hmm searches
    match transfer_hmm_searches(&sqlite) {
        Ok(_) => { },
        Err(e) => panic!("Unable to transfer hmm searches, error: {}", e)
    };
    info!("Successfully transferred SQLite database to RocksDB.  If desired, you may now delete the SQLite database from your hard drive.");

}

fn transfer_hmm_searches(sqlite: &Sqlite) -> Result<(), Error> {

    // Get total est sequences
    let total: u64 = sqlite.get_total("orthograph_hmmsearch");
    info!("Transferring {} hmm searches to RocksDB", total);

    // Execute SQL
    let mut stmt = sqlite.conn.prepare("SELECT s.*,e.header,e.type FROM orthograph_hmmsearch s, orthograph_ests e WHERE s.target = e.digest ORDER BY id").unwrap();
    let mut rows = match stmt.query([]) {
        Ok(r) => r,
        Err(e) => panic!("Unable to execute SQL statement to retrieve est sequencs, error: {}", e)
    };

    // Go through rows
    let mut x = 0;
    loop {

        // Get row
        let row = match rows.next()? {
            Some(r) => r,
            None => break
        };

        // Get blast results
        let hmmsearch_id: u32 = row.get(0)?;
        let blast_results = sqlite.get_blast_results(&hmmsearch_id).unwrap();
        let header: String = row.get(13)?;

        // Set hmm search
        let hmm_search = HmmSearch {
            taxid: row.get(1)?,
            gene: row.get(2)?,
            header: header.replace(" ", "_").to_string(),
            score: row.get(4)?,
            evalue: row.get(5)?,
            log_evalue: row.get(6)?,
            env_start: row.get(7)?,
            env_end: row.get(8)?,
            ali_start: row.get(9)?,
            ali_end: row.get(10)?,
            hmm_start: row.get(11)?,
            hmm_end: row.get(12)?,
            seq_type: row.get(14)?,
            blast: blast_results
        };

        // Set key and json
        let key = format!("hmmsearch:{}", x);
        let json = serde_json::to_string(&hmm_search).unwrap();

        // Add to RocksDB
        ROCKSDB.put(&key, &json);
        x += 1;
        if x % 1000 == 0 {
            print!(".");
            io::stdout().flush().unwrap();
        }
    }

    // Finish
        println!("");
    info!("Successfully transferred total of {} hmm searches", x + 1);

    Ok(())
}

fn transfer_est_sequences(sqlite: &Sqlite) -> Result<(), Error> {

    // Get total est sequences
    let total_est: u64 = sqlite.get_total("orthograph_ests");
    info!("Transferring {} est sequences to RocksDB", total_est);

    // Execute sql statement
    let mut stmt = sqlite.conn.prepare("SELECT id,type,header,sequence FROM orthograph_ests ORDER BY id").unwrap();
    let mut rows = match stmt.query([]) {
        Ok(r) => r,
        Err(e) => panic!("Unable to execute SQL statement to retrieve est sequencs, error: {}", e)
    };

    // Go through rows
    let mut x = 0;
    loop {

        // Get row
        let row = match rows.next()? {
            Some(r) => r,
            None => break
        };
        let header: String = row.get(2)?;

            // Define sequences
        let seq = Sequence {
            seq_type: row.get(1)?,
            sequence: row.get(3)?
        };

        // Add to RocksDB
        ROCKSDB.put(&header.replace(" ", "_"), &seq.sequence);
        if x % 10000 == 0 {
            print!(".");
            io::stdout().flush().unwrap();
        }
        x += 1;
    }
    println!("");
    info!("Successfully saved {} est sequences to RocksDB", x);

    Ok(())
}



