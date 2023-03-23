use clap::{App, Arg};
use std::ops::Range;

pub struct Args {
    pub action: String,
    pub rocksdb: String,
    pub input_dir: String,
    pub sqlite_file: String,
    pub header: String,
    pub seq_type: String,
    pub start: u32,
    pub limit: i32,
    pub coords: Range<usize>,
}

impl Args {
    pub fn new() -> Self {
        // Specify cli arguments
        let matches = App::new("biodb")
            .version("0.1")
            .author("Matt Dizak <matt@apexpl.io>")
            .about("Retrieve the sequences you need.")
            .arg(Arg::with_name("action")
                .short('a')
                .long("action")
                .takes_value(true)
                .help("The action to perform, either 'upgrade-db', 'get-sequence', 'get-hmmsearch' or 'get-hmmsearches'"))
            .arg(Arg::with_name("input")
                .short('i')
                .long("input")
                .takes_value(true)
                .help("Input directory (orthograph_results/species_name)"))
            .arg(Arg::with_name("header")
                .short('h')
                .long("header")
                .takes_value(true)
                .help("Base header of sequence to retrieve."))
            .arg(Arg::with_name("type")
                .short('t')
                .long("type")
                .takes_value(true)
                .help("Type of sequence to retrieve (nt, aa)"))
            .arg(Arg::with_name("coords")
                .short('c')
                .long("coords")
                .takes_value(true)
                .help("Optional coordinates to retrieve, fomatted as START-END (eg. 5-61)"))
            .arg(Arg::with_name("start")
                .short('s')
                .long("start")
                .takes_value(true)
                .help("For action 'get-hmmsearches', and where in the result set to start.  Defaults to 0."))
            .arg(Arg::with_name("limit")
                .short('l')
                .long("limit")
                .takes_value(true)
                .help("For 'get-hmmsearches' action, and the number of hmm searches to return."))
            .get_matches();

        // Get args
        let action = matches.value_of("action").unwrap_or("");
        let input_dir = matches.value_of("input").unwrap_or("/home/boxer/devel/clients/kevin/test-data/old/Syrphidae/orthograph_results/210222/AB34950737_S20_R1.fa");
        let header = matches.value_of("header").unwrap_or("");
        let seq_type = matches.value_of("type").unwrap_or("nt");
        let start = matches.value_of("start").unwrap_or("0");
        let limit = matches.value_of("limit").unwrap_or("500");

        // Get coords
        let coords_str = matches.value_of("coords").unwrap_or("0-0").split("-");
        let coords = coords_str
            .into_iter()
            .map(|c| c.parse::<usize>().unwrap())
            .collect::<Vec<usize>>();

        // Parse input dir
        let parts = input_dir.split("/");
        let sqlite_file = format!(
            "{}/{}.sqlite",
            input_dir,
            parts
                .last()
                .unwrap()
                .trim_end_matches("/")
                .trim_end_matches(".fa")
                .to_string()
        );

        // Return
        Self {
            action: action.to_string(),
            input_dir: input_dir.trim_end_matches("/").to_string(),
            rocksdb: format!("{}/rocksdb", input_dir.trim_end_matches("/").to_string()).to_string(),
            sqlite_file: sqlite_file,
            header: header.to_string(),
            seq_type: seq_type.to_string(),
            start: start.parse::<u32>().unwrap(),
            limit: limit.parse::<i32>().unwrap(),
            coords: coords[0]..coords[1],
        }
    }
}
