#![allow(warnings)]
use crate::reporter::Reporter;
use env_logger::{Builder, Target};
use log::{info, LevelFilter};
use std::io::Write;
use std::time::Instant;
use biotools::CONFIG;

mod reporter;
mod algorithms;
mod output;
mod stats;
mod models;
mod temp_tables;

fn main() {

    // Greeting
    greeting();

    // Initialize logger
    init_logger();
    let start_time = Instant::now();

    // Process reporter
    let reporter = Reporter::new();
    reporter.process();

    // Give processing time
    let elapsed = start_time.elapsed();
    info!("Completed processing in {:?} seconds.", elapsed.as_secs());
}

fn greeting() {

    println!("Orthograph: Orthology prediction using a Graph-based,");
    println!("Reciprocal Approach with Profile Hidden Markov models");
    println!("      Originally in Perl by Malte Petersen <mptrsen@uni-bonn.de> (2015)");
    println!("      Converted to Rust by Matt Dizak <matt@apexpl.io> (May 2022)");
    println!("      Version: {}", env!("CARGO_PKG_VERSION"));
    println!("");
}

fn init_logger() {

    // Get log level
    let mut log_level = LevelFilter::Warn;
    if CONFIG.log.verbose == true {
        log_level = LevelFilter::Debug;
    } else if !CONFIG.log.quiet {
        log_level = LevelFilter::Info;
    }

    // Init logger
    Builder::new()
        .format(|buf, record| {
            writeln!(buf, "{}: {}", record.level(), record.args())
        }).filter(None, log_level).target(Target::Stdout).init();

    info!("Initialized logging at {}", CONFIG.log.logfile);
    }



