extern crate serde;

use crate::database::Database;
use crate::{BIODB_ARGS, ROCKSDB};
use std::io::{self, Write};
use std::string::String;
use log::error;

pub fn get() {

    // Get hmm search
    let key = format!("hmmsearch:{}", BIODB_ARGS.header);
    let json: String = match ROCKSDB.get(&key) {
        Some(r) => r,
        None => {
            error!("No hmm search exists with the id# {}", BIODB_ARGS.header);
            return;
        }
    };

    // Print result
    io::stdout().write_all(json.as_bytes()).unwrap();
}

pub fn get_multi() {

    // Initialize
    let mut iter = ROCKSDB.db.raw_iterator();
    let key = format!("hmmsearch:{}", BIODB_ARGS.start);
    let mut x: i32 = 0;
    io::stdout().write_all("[".as_bytes()).unwrap();

    // Go through rows
    iter.seek(&key.as_bytes());
    while iter.valid() {

        let json = match iter.value() {
            Some(r) => String::from_utf8(r.to_vec()).unwrap(),
            None => break
        };

        if x > 0 {
            io::stdout().write_all(",".as_bytes()).unwrap();
        }
        io::stdout().write_all(json.as_bytes()).unwrap();
        io::stdout().flush().unwrap();

        if BIODB_ARGS.limit > -1 && x >= BIODB_ARGS.limit {
            break;
        }
        iter.next();
        x += 1;
    }

    io::stdout().write_all("]".as_bytes()).unwrap();
}


