use rusqlite::Connection;
use std::collections::HashMap;
use biotools::CONFIG;

pub fn setup(conn: &Connection) {

    // Ensure we're down
    teardown(&conn);
    let mut table_sql = HashMap::new();

    // Hits table
    table_sql.insert("hits", format!("CREATE TABLE {}_hits (
        id INTEGER PRIMARY KEY,
        is_overlap BOOLEAN NOT NULL DEFAULT true,
        hmmsearch_id unsigned integer not null,
        taxid unsigned integer not null,
        aaseq_id unsigned integer not null,
        ntseq_id unsigned integer not null,
        blast_target unsigned integer not null,
        gene_id VARCHAR(255) NOT NULL,
        score DOUBLE NOT NULL,
        digest VARCHAR(32) NOT NULL,
        evalue VARCHAR(8) NOT NULL,
        hmm_start UNSIGNED INTEGER NOT NULL,
        hmm_end UNSIGNED INTEGER NOT NULL,
        ali_start UNSIGNED INTEGER NOT NULL,
        ali_end UNSIGNED INTEGER NOT NULL,
        env_start UNSIGNED INTEGER NOT NULL,
        env_end UNSIGNED INTEGER NOT NULL,
        blast_start UNSIGNED INTEGER NOT NULL,
        blast_end UNSIGNED INTEGER NOT NULL,
        header_base VARCHAR(255) NOT NULL,
        header_full VARCHAR(255) NOT NULL,
        header_revcomp BOOLEAN NOT NULL DEFAULT false,
        header_translate UNSIGNED INTEGER NOT NULL,
        non_orf_sequence BLOB NOT NULL
    )", CONFIG.db.table_prefix));

    // orf table
    table_sql.insert("orf", format!("CREATE TABLE {}_orf (
        id INTEGER PRIMARY KEY,
        hit_id UNSIGNED INTEGER NOT NULL,
        taxid UNSIGNED INTEGER NOT NULL,
        cdna_start UNSIGNED INTEGER NOT NULL,
        cdna_end UNSIGNED INTEGER NOT NULL,
        aa_start UNSIGNED INTEGER NOT NULL,
        aa_end UNSIGNED INTEGER NOT NULL,
        cdna_start_transcript UNSIGNED INTEGER NOT NULL,
        cdna_end_transcript UNSIGNED INTEGER NOT NULL,
        aa_start_transcript UNSIGNED INTEGER NOT NULL,
        aa_end_transcript UNSIGNED INTEGER NOT NULL,
        aa_start_hmm UNSIGNED INTEGER NOT NULL,
        aa_end_hmm UNSIGNED INTEGER NOT NULL,
        translated_seq BLOB NOT NULL,
        cdna_seq MEDIUM BLOB NOT NULL
    )", CONFIG.db.table_prefix));

    // Go through and create tables
    for table_name in table_sql.keys() {
        match conn.execute(&table_sql[table_name], []) {
            Ok(res) => res,
            Err(e) => panic!("Unable to create temporary database table '{}' with error: {}", table_name, e)
        };
    }

}

pub fn teardown(conn: &Connection) {

    // Drop tables
    let tables = vec!["hits", "orf"];
    for table in tables {
        let sql = format!("DROP TABLE IF EXISTS {}_{}", CONFIG.db.table_prefix, table);
        match conn.execute(&sql, []) {
            Ok(res) => res,
            Err(e) => panic!("Unable to drop temporary database table {} with error: {}", table, e)
        };
    }

}


