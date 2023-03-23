use crate::reporter::ReporterKit;
use biotools::db::sqlite::{TBL_HITS, TBL_LOGS, TBL_SEQUENCE_PAIRS, TBL_TAXA};
use biotools::CONFIG;
use log::info;
use rusqlite::{Error, Statement};
use std::fs::File;
use std::io::Write;
use std::path::Path;

struct GeneCount {
    gene_id: String,
    count: u32,
}

struct CoreSequence {
    gene_id: String,
    taxa_name: String,
    header: String,
    sequence: String,
}

struct Sequence {
    gene_id: String,
    header: String,
    is_revcomp: u8,
    translate: u8,
    taxa_name: String,
    aa_start: u16,
    aa_end: u16,
    cdna_start: u16,
    cdna_end: u16,
    aa_seq: String,
    cdna_seq: String,
}

pub fn run(kit: &ReporterKit) -> Result<bool, Error> {
    // Prepare sql
    let sql = format!(
        "SELECT gene_id,count(gene_id) FROM {} GROUP BY gene_id ORDER BY gene_id",
        *TBL_HITS
    );
    let mut stmt = match kit.memdb.prepare(&sql) {
        Ok(r) => r,
        Err(e) => panic!(
            "Unable to prepare sql to select genes for finalization, error: {}",
            e
        ),
    };
    info!("Writing sequence files of all genes");

    // GO through rows
    let mut rows = match stmt.query([]) {
        Ok(r) => r,
        Err(e) => panic!(
            "Unable to execute sql to gather genes during writing of sequence files, error: {}",
            e
        ),
    };
    // GO through rows
    loop {
        // Get row
        let row = match rows.next()? {
            Some(r) => r,
            None => break,
        };

        // Set variables
        let gene_id: String = row.get(0)?;
        let count: u32 = row.get(1)?;

        // Write sequences files, if needed
        if count > 0 {
            save_gene(&kit, &gene_id);
        }
    }

    Ok(true)
}

fn save_gene(kit: &ReporterKit, gene_id: &String) {
    // Prepare
    let tmp_gene = format!("{}", gene_id);
    let (mut aa_fh, mut nt_fh) = create_files(tmp_gene);

    // Save core sequences
    write_core_sequences(&kit, &gene_id, &mut aa_fh, "aa".to_string());
    write_core_sequences(&kit, &gene_id, &mut nt_fh, "nt".to_string());

    // Save sequences
    match write_sequences(&kit, &gene_id, &mut aa_fh, "aa".to_string()) {
        Ok(r) => r,
        Err(e) => panic!("Unable to save sequences, {}", e),
    };
    write_sequences(&kit, &gene_id, &mut nt_fh, "nt".to_string());
}

fn write_core_sequences(
    kit: &ReporterKit,
    gene_id: &String,
    mut fh: &mut File,
    seq_type: String,
) -> Result<bool, Error> {
    // Execute sql
    let mut stmt = prepare_core_sql(&kit, seq_type);
    let mut rows = match stmt.query([&gene_id]) {
        Ok(r) => r,
        Err(e) => panic!(
            "Unable to execute sql to select hits while writing final sequence files, error: {}",
            e
        ),
    };

    // GO through rows
    loop {
        // Get row
        let row = match rows.next()? {
            Some(r) => r,
            None => break,
        };

        // Define sequence
        let seq = CoreSequence {
            gene_id: row.get(0)?,
            taxa_name: row.get(1)?,
            header: row.get(2)?,
            sequence: row.get(3)?,
        };

        // Format header
        let header = vec![
            seq.gene_id,
            seq.taxa_name,
            seq.header,
            format!("1-{}", seq.sequence.len().to_string()).to_string(),
            ".".to_string(),
            ".".to_string(),
        ]
        .join(&CONFIG.search.header_seperator);

        // Save to aa file
        let line = format!(">{}\n{}\n", header, seq.sequence);
        fh.write_all(&line.as_bytes())
            .expect("Unable to write to aa results file");
    }

    Ok(true)
}

fn write_sequences(
    kit: &ReporterKit,
    gene_id: &String,
    mut fh: &mut File,
    seq_type: String,
) -> Result<bool, Error> {
    // Execute sql
    let mut stmt = prepare_sequence_sql(&kit);
    let mut rows = match stmt.query([&gene_id]) {
        Ok(r) => r,
        Err(e) => panic!(
            "Unable to execute sql to select hits while writing final sequence files, error: {}",
            e
        ),
    };

    // GO through rows
    loop {
        // Get row
        let row = match rows.next()? {
            Some(r) => r,
            None => break,
        };

        // Define sequence
        let seq = Sequence {
            gene_id: row.get(0)?,
            header: row.get(1)?,
            is_revcomp: row.get(2)?,
            translate: row.get(3)?,
            taxa_name: row.get(4)?,
            aa_start: row.get(5)?,
            aa_end: row.get(6)?,
            cdna_start: row.get(7)?,
            cdna_end: row.get(8)?,
            aa_seq: row.get(9)?,
            cdna_seq: row.get(10)?,
        };

        // Get rf
        let rf = get_rf(&seq.is_revcomp, &seq.translate, &seq_type);

        // Get header
        let header = vec![
            seq.gene_id,
            format!("{}", CONFIG.report.species_name),
            seq.header,
            format!("{}-{}", seq.aa_start, seq.aa_end),
            rf.to_string(),
            seq.taxa_name,
        ]
        .join(&CONFIG.search.header_seperator);

        // Get sequence
        let sequence = if seq_type == "nt".to_string() {
            seq.cdna_seq
        } else {
            seq.aa_seq
        };

        // Save to aa file
        let line = format!(">{}\n{}\n", header, sequence);
        fh.write_all(&line.as_bytes())
            .expect("Unable to write to aa results file");
    }

    Ok(true)
}

fn create_files(gene_id: String) -> (File, File) {
    // Open aa file
    let aa_filename = format!("{}/aa/{}.aa.fa", CONFIG.report.output_dir, gene_id);
    let aa_path = Path::new(&aa_filename);
    let mut aa_fh = match File::create(&aa_path) {
        Ok(res) => res,
        Err(e) => panic!(
            "Unable to open file for writing, {}, error: {}",
            aa_filename, e
        ),
    };

    // Open nt file
    let nt_filename = format!("{}/nt/{}.nt.fa", CONFIG.report.output_dir, gene_id);
    let nt_path = Path::new(&nt_filename);
    let mut nt_fh = match File::create(&nt_path) {
        Ok(res) => res,
        Err(e) => panic!(
            "Unable to open file for writing, {}, error: {}",
            nt_filename, e
        ),
    };

    // Return
    (aa_fh, nt_fh)
}

fn prepare_core_sql(kit: &ReporterKit, seq_type: String) -> Statement {
    // Initialize
    let seq_table = format!("input.{}_{}seqs", CONFIG.db.table_prefix, seq_type);

    // Define sql
    let sql = format!(
        "SELECT l.ortholog_gene_id, t.name, a.header, a.sequence 
        FROM {} a, {} t, {} p, {} l 
    WHERE 
        l.ortholog_gene_id = ? AND  
        p.{}_seq = a.id AND 
        l.sequence_pair = p.id AND 
        a.taxid = t.id 
        ORDER BY t.name, a.header
    ",
        seq_table, *TBL_TAXA, *TBL_SEQUENCE_PAIRS, *TBL_LOGS, seq_type
    );

    // Prepare sql
    let stmt = match kit.memdb.prepare(&sql) {
        Ok(r) => r,
        Err(e) => panic!("Unable to prepare sql statement to retrieve hits to save final sequence files, error: {}", e)
    };

    stmt
}

fn prepare_sequence_sql(kit: &ReporterKit) -> Statement {
    // SDefine sql
    let sql = format!(
        "SELECT 
        h.gene_id,
        h.header_base,
        h.header_revcomp,
        h.header_translate,
        t.name,
        o.aa_start_transcript,
        o.aa_end_transcript,
            o.cdna_start_transcript,
            o.cdna_end_transcript,
            o.translated_seq,
            o.cdna_seq 
        FROM 
            {} h, {}_orf o, {} t  
        WHERE
            h.gene_id = ? AND 
            h.id = o.hit_id AND 
            t.id = o.taxid AND 
            LENGTH(o.translated_seq) >= {} AND 
            LENGTH(o.cdna_seq) >= {} 
            GROUP BY o.translated_seq ORDER BY h.id 
        ",
        *TBL_HITS,
        CONFIG.db.table_prefix,
        *TBL_TAXA,
        CONFIG.search.min_transcript_length,
        CONFIG.search.min_transcript_length
    );

    // Prepare sql
    let stmt = match kit.memdb.prepare(&sql) {
        Ok(r) => r,
        Err(e) => panic!(
            "Unable to prepare SQL to retrieve hits during sequence finalization, error: {}",
            e
        ),
    };

    stmt
}

fn get_rf(revcomp: &u8, translate: &u8, seq_type: &str) -> String {
    // Check for nt
    if seq_type == "nt" {
        return String::from(".");
    }

    // Get rf
    if revcomp == &1 && translate > &0 {
        return format!("[revcomp]:[translate({})]", translate).to_string();
    } else if revcomp == &1 {
        return String::from("[revcomp]");
    } else if translate > &0 {
        return format!("[translate({})]", translate).to_string();
    }

    "".to_string()
}
