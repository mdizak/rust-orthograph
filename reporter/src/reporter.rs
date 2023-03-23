use crate::algorithms::{
    env_pseudo_master, extract_reciprocal_hits, frameshift_correction, hmm_overlap,
};
use crate::output::{save_brh_files, save_sequence_files};
use crate::temp_tables;
use biotools::db::sqlite::Sqlite;
use biotools::CONFIG;
use log::info;
use rusqlite::Connection;
use std::collections::HashMap;

pub struct Reporter {}

pub struct ReporterKit {
    pub db: Sqlite,
    pub memdb: Connection,
    pub species_id: u32,
    pub set_id: u32,
    pub reference_taxa: Vec<String>,
    pub aaseq_by_gene: HashMap<String, Vec<u32>>,
}

impl Reporter {
    pub fn new() -> Self {
        Self {}
    }

    pub fn process(self) {
        // Prepare environment
        self.prepare();

        // Initialize
        let kit = self.initialize();

        // Extract reciprocal hits
        let mut stats = extract_reciprocal_hits::run(&kit)
            .expect("Error occured while trying to extract reciprocal hits.");

        // Env pseudo master checks
        env_pseudo_master::check(&kit, &mut stats);

        // Hmm overlap checks
        hmm_overlap::check(&kit, &mut stats);

        // Save brh files
        save_brh_files::save(&kit, &mut stats);

        // Frameshift correction
        frameshift_correction::run(&kit, &mut stats);

        // Save sequence files
        save_sequence_files::run(&kit);

        // Write report
        stats.write_report();
    }

    fn prepare(&self) {
        // Check output directory, and create if necessary
        biotools::io::create_dir(&CONFIG.report.output_dir);

        // Create log directory, if needed
        let logdir = format!("{}/log", CONFIG.report.output_dir);
        biotools::io::create_dir(&logdir);

        // Clear and re-create /aa/ and /nt/ directories, if needed
        if CONFIG.switch.clear_files == true {
            let aa_dir = format!("{}/aa", &CONFIG.report.output_dir);
            let nt_dir = format!("{}/nt", &CONFIG.report.output_dir);
            biotools::io::recreate_dir(&aa_dir);
            biotools::io::recreate_dir(&nt_dir);
        }
    }

    fn initialize(&self) -> ReporterKit {
        // Connect to SQLite, get species id
        let db = Sqlite::new();
        let species_id: u32 = match db.get_species_id() {
            Ok(res) => res,
            Err(e) => panic!(
                "Unable to determine id# for species {}, error: {}",
                &CONFIG.report.species_name, e
            ),
        };
        info!(
            "Got species id# {} for species name {}",
            species_id, CONFIG.report.species_name
        );

        // Get set id
        let set_id: u32 = match db.get_set_id() {
            Ok(res) => res,
            Err(e) => panic!(
                "Unable to determine id# for set {}, error: {}",
                &CONFIG.report.set_name, e
            ),
        };
        info!(
            "Got set id# {} for set name {}",
            set_id, CONFIG.report.set_name
        );

        // Get referenced taxa
        let ref_taxa: Vec<String> = if CONFIG.report.reference_taxa.len() > 0 {
            CONFIG
                .report
                .reference_taxa
                .split(",")
                .map(str::to_string)
                .filter(|i| i != "")
                .collect::<Vec<String>>()
        } else {
            db.get_reference_taxa(&set_id).unwrap()
        };
        info!(
            "Obtained {} reference taxa to use.",
            ref_taxa.len().to_string()
        );

        // Get aa sequences in set
        let aaseq = match db.get_aaseq_in_set(&set_id) {
            Ok(res) => res,
            Err(e) => panic!("Unable to retrieve aa sequences within set, error: {}", e),
        };
        info!(
            "Obtained total of {} genes with aa sequences for reporting.",
            aaseq.len().to_string()
        );

        // Return
        ReporterKit {
            db: db,
            memdb: self.open_memdb(),
            species_id: species_id,
            set_id: set_id,
            reference_taxa: ref_taxa,
            aaseq_by_gene: aaseq,
        }
    }

    fn open_memdb(&self) -> Connection {
        // Connect to in-memory SQLite database
        let memdb = Connection::open_in_memory()
            //let memdb = Connection::open("/home/boxer/devel/clients/kevin/rust/report.db")
            .expect("Unable to open in-memory SQLite database.");

        // Attach database
        let sql = format!("ATTACH '{}' AS out", CONFIG.db.reporter_sqlite_file);
        memdb
            .execute(&sql, [])
            .expect("Unable to attach output SQLite database to in-memory database");

        // Attach input database
        let sql = format!("ATTACH '{}' AS input", CONFIG.db.sqlite_file);
        memdb
            .execute(&sql, [])
            .expect("Unable to attach input SQLite database to in-memory database");

        // Create temporary database tables
        temp_tables::setup(&memdb);
        info!("Created temporary database tables...");

        // Return
        memdb
    }
}
