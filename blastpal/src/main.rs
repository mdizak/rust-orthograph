use biotools::settings::Settings;
use crate::blastpal::Blastpal;
use lazy_static::lazy_static;
use env_logger::{Builder, Target};
use log::{info, LevelFilter};
use std::io::Write;
use std::time::Instant;

mod blastpal;
mod blast;
mod models;

lazy_static! {
    pub static ref CONFIG: Settings = Settings::new();
}


fn main() {

    // Init logger
    init_logger();
    let start_time = Instant::now();

    // Process
    let blastpal = Blastpal::new();
    match blastpal.process() {
        Ok(res) => res,
        Err(e) => panic!("An error occured while processing: {}", e)
    };

    // Give processing time
    let elapsed = start_time.elapsed();
    info!("Completed processing in {:?} seconds.", elapsed.as_secs());
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



