use chrono::prelude::*;
use config::Config;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use std::string::String;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Database {
    pub sqlite_file: String,
    pub reporter_sqlite_file: String,
    pub table_prefix: String,
}

pub struct Search {
    pub hmmsearch_threshold: u16,
    pub blast_threshold: u16,
    pub hmmsearch_evalue_threshold: f32,
    pub blast_evalue_threshold: f32,
    pub env_overlap_threshold: f32,
    pub env_score_discard_threshold: f32,
    pub hmm_overlap_threshold: f32,
    pub hmm_score_discard_threshold: f32,
    pub max_blast_searches: u16,
    pub max_blast_hits: u16,
    pub num_threads: u8,
    pub min_transcript_length: u16,
    pub min_overlap: f32,
    pub fill_with_x: bool,
    pub substitute_u_with: String,
    pub header_seperator: String,
    pub max_mismatches: u16,
}

pub struct Switch {
    pub brh_only: bool,
    pub frameshift_correction: bool,
    pub extend_orf: bool,
    pub strict_search: bool,
    pub clear_database: bool,
    pub clear_files: bool,
    pub enable_env_overlap: bool,
    pub enable_hmm_overlap: bool,
}

pub struct Log {
    pub verbose: bool,
    pub quiet: bool,
    pub logfile: String,
}

pub struct Report {
    pub input_file: String,
    pub blastdb: String,
    pub species_name: String,
    pub set_name: String,
    pub sets_dir: String,
    pub output_dir: String,
    pub reference_taxa: String,
    pub cog_list_file: String,
    pub wanted_genes: Vec<String>,
}

pub struct Programs {
    pub sqlite: String,
    pub alignment: String,
    pub hmmbuild: String,
    pub makeblastdb: String,
    pub translate: String,
    pub hmmsearch: String,
    pub blast: String,
    pub exonerate: String,
}

pub struct Settings {
    pub db: Database,
    pub search: Search,
    pub switch: Switch,
    pub log: Log,
    pub report: Report,
    pub programs: Programs,
}

impl Settings {
    pub fn new() -> Self {
        // Load config file
        let config = Settings::load_config_file();

        // Database
        let database = Database {
            sqlite_file: Settings::get_var(&config, &"sqlite-database"),
            reporter_sqlite_file: format!(
                "{}/{}.sqlite",
                config["output-directory"], config["species-name"]
            ),
            table_prefix: Settings::get_var(&config, "dbtable-prefix"),
        };

        // Search
        let search = Search {
            hmmsearch_threshold: Settings::get_var(&config, "hmmsearch-score-threshold")
                .parse::<u16>()
                .unwrap(),
            blast_threshold: Settings::get_var(&config, "blast-score-threshold")
                .parse::<u16>()
                .unwrap(),
            hmmsearch_evalue_threshold: Settings::get_var(&config, "hmmsearch-evalue-threshold")
                .parse::<f32>()
                .unwrap(),
            blast_evalue_threshold: Settings::get_var(&config, "blast-evalue-threshold")
                .parse::<f32>()
                .unwrap(),
            env_overlap_threshold: Settings::get_var(&config, "env-overlap-threshold")
                .parse::<f32>()
                .unwrap(),
            env_score_discard_threshold: Settings::get_var(&config, "env-score-discard-threshold")
                .parse::<f32>()
                .unwrap(),
            hmm_overlap_threshold: Settings::get_var(&config, "hmm-overlap-threshold")
                .parse::<f32>()
                .unwrap(),
            hmm_score_discard_threshold: Settings::get_var(&config, "hmm-score-discard-threshold")
                .parse::<f32>()
                .unwrap(),
            max_blast_searches: Settings::get_var(&config, "max-blast-searches")
                .parse::<u16>()
                .unwrap(),
            max_blast_hits: Settings::get_var(&config, "max-blast-hits")
                .parse::<u16>()
                .unwrap(),
            num_threads: Settings::get_var(&config, "num-threads")
                .parse::<u8>()
                .unwrap(),
            min_transcript_length: Settings::get_var(&config, "minimum-transcript-length")
                .parse::<u16>()
                .unwrap(),
            min_overlap: Settings::get_var(&config, "orf-overlap-minimum")
                .parse::<f32>()
                .unwrap(),
            fill_with_x: Settings::get_var(&config, "fill-with-x")
                .parse::<bool>()
                .unwrap(),
            substitute_u_with: Settings::get_var(&config, "substitute-u-with"),
            header_seperator: Settings::get_var(&config, "header-separator"),
            max_mismatches: Settings::get_var(&config, "max-reciprocal-mismatches")
                .parse::<u16>()
                .unwrap(),
        };

        // Switch
        let switch = Switch {
            brh_only: Settings::get_bool(&config, "brh-only"),
            frameshift_correction: Settings::get_bool(&config, "frameshift-correction"),
            extend_orf: Settings::get_bool(&config, "extend-orf"),
            strict_search: Settings::get_bool(&config, "strict-search"),
            clear_database: Settings::get_bool(&config, "clear-database"),
            clear_files: Settings::get_bool(&config, "clear-files"),
            enable_env_overlap: Settings::get_bool(&config, "enable-env-overlap"),
            enable_hmm_overlap: Settings::get_bool(&config, "enable-hmm-overlap"),
        };

        // Log
        let log = Log {
            verbose: Settings::get_bool(&config, "verbose"),
            quiet: Settings::get_var(&config, "quiet").parse::<bool>().unwrap(),
            logfile: Settings::get_var(&config, "logfile"),
        };

        // Report
        let report = Report {
            input_file: Settings::get_var(&config, "input-file"),
            blastdb: Settings::get_var(&config, "blastdb"),
            species_name: Settings::get_var(&config, "species-name"),
            set_name: Settings::get_var(&config, "ortholog-set"),
            output_dir: Settings::get_var(&config, "output-directory")
                .trim_end_matches("/")
                .to_string(),
            reference_taxa: Settings::get_var(&config, "reference-taxa"),
            sets_dir: Settings::get_var(&config, "sets-dir")
                .trim_end_matches("/")
                .to_string(),
            cog_list_file: Settings::get_var(&config, "cog-list-file"),
            wanted_genes: Settings::get_wanted_genes(&Settings::get_var(&config, "cog-list-file")),
        };

        // Programs
        let programs = Programs {
            sqlite: Settings::get_var(&config, "sqlite-program"),
            alignment: Settings::get_var(&config, "alignment-program"),
            hmmbuild: Settings::get_var(&config, "hmmbuild-program"),
            makeblastdb: Settings::get_var(&config, "makeblastdb-program"),
            translate: Settings::get_var(&config, "translate-program"),
            hmmsearch: Settings::get_var(&config, "hmmsearch-program"),
            blast: Settings::get_var(&config, "blast-program"),
            exonerate: Settings::get_var(&config, "exonerate-program"),
        };

        // Return
        Self {
            db: database,
            search: search,
            switch: switch,
            log: log,
            report: report,
            programs: programs,
        }
    }

    pub fn get_wanted_genes(cogfile: &String) -> Vec<String> {
        // Check file exists
        let mut genes: Vec<String> = Vec::new();
        if !Path::new(cogfile).exists() {
            return genes;
        }

        // Open file
        let fh = match File::open(&cogfile) {
            Ok(res) => res,
            Err(e) => panic!("Unable to open cog-list-file at {}, error: {}", cogfile, e),
        };
        let lines = io::BufReader::new(fh).lines();

        // Go through lines
        for line in lines {
            if let Ok(gene) = line {
                let fgene = gene.trim_end();
                if fgene != "" {
                    genes.push(fgene.to_string());
                }
            }
        }

        // Return
        genes
    }

    /**
     * Get single variable from config hashmap
     */
    fn get_var(config: &HashMap<String, String>, name: &str) -> String {
        let key = String::from(name);
        config.get(&key).unwrap().to_string()
    }

    fn get_bool(config: &HashMap<String, String>, name: &str) -> bool {
        let key = String::from(name);
        let value: &str = config.get(&key).unwrap();
        let res = match value {
            "true" => true,
            "false" => false,
            "1" => true,
            "0" => false,
            _ => false,
        };
        res
    }

    /**
     * Load config file, return hashmap of contents
     */
    fn load_config_file() -> HashMap<String, String> {
        // Get config file
        let config_file = if std::env::args().len() > 2 && std::env::args().nth(1).unwrap() == "-c"
        {
            std::env::args().nth(2).unwrap().to_string()
        } else {
            "config.ini".to_string()
        };

        // Check file exists
        if !Path::new(&config_file).exists() {
            panic!(
                "No {} configuration file exists within this directory.",
                config_file
            );
        }

        // Read file
        let settings = Config::builder()
            .add_source(config::File::with_name(&config_file))
            .build()
            .unwrap()
            .try_deserialize::<HashMap<String, String>>();

        // Get default values
        let mut config = Settings::get_default_values();

        // Go through settings, add to config hashmap
        for parent in &settings {
            for (k, v) in parent {
                config.insert(k.to_string(), v.to_string());
            }
        }

        // Validate and prepare
        Settings::validate(&config);

        // Return
        config
    }

    /**
     * Get default values
     */
    fn get_default_values() -> HashMap<String, String> {
        // Get default log file
        let lt = Local::now();
        let logfile: String = format!(
            "orthograph-reporter-{}-{:02}-{:02}_{:02}:{:02}.log",
            lt.year(),
            lt.month(),
            lt.day(),
            lt.hour(),
            lt.minute()
        );

        let config = HashMap::from([
            (String::from("blastdb"), String::from("")),
            (String::from("dbtable-prefix"), String::from("orthograph")),
            (String::from("sqlite-program"), String::from("sqlite3")),
            (String::from("ortholog-set"), String::from("test_set")),
            (String::from("hmmbuild-program"), String::from("hmmbuild")),
            (
                String::from("makeblastdb-program"),
                String::from("makeblastdb"),
            ),
            (
                String::from("translate-program"),
                String::from("fastatranslate"),
            ),
            (String::from("hmmsearch-program"), String::from("hmmsearch")),
            (String::from("blast-program"), String::from("blast")),
            (String::from("exonerate-program"), String::from("exonerate")),
            (
                String::from("alignment-program"),
                String::from("mafft-linsi"),
            ),
            (
                String::from("hmmsearch-score-threshold"),
                String::from("10"),
            ),
            (String::from("blast-score-threshold"), String::from("10")),
            (String::from("env-overlap-threshold"), String::from("20")),
            (
                String::from("env-score-discard-threshold"),
                String::from("70"),
            ),
            (String::from("hmm-overlap-threshold"), String::from("50")),
            (
                String::from("hmm-score-discard-threshold"),
                String::from("1.5"),
            ),
            (
                String::from("hmmsearch-evalue-threshold"),
                String::from("1e-05"),
            ),
            (
                String::from("blast-evalue-threshold"),
                String::from("1e-05"),
            ),
            (String::from("max-blast-searches"), String::from("100")),
            (String::from("max-blast-hits"), String::from("100")),
            (String::from("num-threads"), String::from("1")),
            (
                String::from("minimum-transcript-length"),
                String::from("30"),
            ),
            (String::from("fill-with-x"), String::from("false")),
            (String::from("brh-only"), String::from("false")),
            (String::from("frameshift-correction"), String::from("true")),
            (String::from("orf-overlap-minimum"), String::from("0.5")),
            (String::from("strict-search"), String::from("false")),
            (String::from("clear-database"), String::from("true")),
            (String::from("clear-files"), String::from("false")),
            (String::from("enable-env-overlap"), String::from("true")),
            (String::from("enable-hmm-overlap"), String::from("true")),
            (String::from("verbose"), String::from("false")),
            (String::from("quiet"), String::from("false")),
            (String::from("extend-orf"), String::from("false")),
            (String::from("substitute-u-with"), String::from("X")),
            (String::from("header-separator"), String::from("|")),
            (String::from("logfile"), logfile),
            (String::from("reference-taxa"), String::from("")),
            (String::from("cog-list-file"), String::from("")),
            (String::from("max-reciprocal-mismatches"), String::from("0")),
        ]);

        // return
        config
    }

    /**
     * Validate and prepare the config hashmap for struct creation.
     */
    fn validate(config: &HashMap<String, String>) {
        // Check required
        let required = ["output-directory", "species-name"];
        for req in required {
            if !config.contains_key(req) {
                panic!(
                    "No '{}' setting defined within config.ini, which is a required setting.",
                    req
                );
            }
        }

        // Ensure input file exists
        if !Path::new(&config["input-file"]).exists() {
            panic!(
                "The input file does not exist at {}.  Please check the config.ini file.",
                config["input-file"]
            );
        }
    }
}
