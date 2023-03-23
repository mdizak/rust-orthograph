#![allow(warnings)]
use crate::args::Args;
use crate::database::{Database, RocksDB};
use env_logger::{Builder, Target};
use log::{error, LevelFilter};
use std::io::Write;
use lazy_static::lazy_static;

mod args;
mod models;
mod database;
mod sqlite;
mod upgrade_db;
mod sequence;
mod translate;
mod hmmsearch;

lazy_static! {
    pub static ref BIODB_ARGS: Args = Args::new();
    pub static ref ROCKSDB: RocksDB = Database::new();
}

fn main() {

    // Initialize logger
    init_logger();

    // Perform action
    match BIODB_ARGS.action.as_str() {
        "upgrade-db" => upgrade_db::upgrade(),
        "get-sequence" => sequence::get(),
        "get-hmmsearch" => hmmsearch::get(),
        "get-hmmsearches" => hmmsearch::get_multi(),
        _ => error!("Usage: biodb -a (upgrade-db|get-sequence|get-hmmsearch|get-hmmsearches) [-i <INPUT_DIR>] [-h <HEADER|HMM_SEARCH_ID>] [-c S-E] [-t <aa|nt>] [-s <START>] [-l <LIMIT>]]")
    };

}

fn init_logger() {

    // Get log level
    let log_level = LevelFilter::Debug;

    // Init logger
    Builder::new()
        .format(|buf, record| {
            writeln!(buf, "{}: {}", record.level(), record.args())
        }).filter(None, log_level).target(Target::Stdout).init();

}





