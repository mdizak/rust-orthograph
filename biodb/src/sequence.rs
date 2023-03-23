extern crate serde;

use crate::database::Database;
use crate::models::Sequence;
use crate::{BIODB_ARGS, ROCKSDB};
use log::error;
use std::io::{self, Write};
use std::string::String;

pub fn get() {
    // Get sequence
    let mut seq: String = match ROCKSDB.get(&BIODB_ARGS.header) {
        Some(r) => r,
        None => {
            error!("No sequence exists with header: {}", BIODB_ARGS.header);
            return;
        }
    };

    // Extract sequence from json
    //let mut seq: Sequence = serde_json::from_str(&json).unwrap();

    // Check for coords
    if BIODB_ARGS.coords.end > 0 {
        seq = seq[BIODB_ARGS.coords.start..BIODB_ARGS.coords.end].to_string();
    }

    // Translate sequence, if needed
    if BIODB_ARGS.seq_type == "aa".to_string() {
        seq = crate::translate::translate(&seq);
    }

    // Output sequence
    io::stdout().write_all(seq.as_bytes()).unwrap();
}
