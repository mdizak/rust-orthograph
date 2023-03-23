use crate::settings::Settings;
use lazy_static::lazy_static;
use regex::Regex;
use std::collections::HashMap;
use std::ops::Range;

pub mod db;
pub mod io;
pub mod settings;

lazy_static! {
    pub static ref CONFIG: Settings = Settings::new();
    static ref HDR_REGEX: Regex = Regex::new(r"\[(.+?)(\((\d)\))?\]").unwrap();
    static ref HDR_CLEAN_REGEX: Regex = Regex::new(r"\[.*$").unwrap();
}

pub struct FastaResult {
    pub coord_start: u16,
    pub coord_end: u16,
    pub sequence: String,
}

pub fn translate_header(header: &String) -> (String, bool, u8) {
    // Check for extra info
    let (mut is_revcomp, mut translate): (bool, u8) = (false, 0);
    for cap in HDR_REGEX.captures_iter(&header) {
        if cap[1].to_string() == "revcomp" {
            is_revcomp = true;
        }
        if cap[1].to_string() == "translate" {
            translate = cap[3].parse::<u8>().unwrap();
        }
    }

    // Clean right side
    let base_header = HDR_CLEAN_REGEX.replace_all(&header, "");

    // Return
    (
        base_header.trim_end_matches(":").trim_end().to_string(),
        is_revcomp,
        translate,
    )
}

pub fn format_header(header: &String, revcomp: &u8, translate: &u8) -> String {
    // Format as necessary
    if revcomp == &1 && translate > &0 {
        return format!("{} [revcomp]:[translate({})]", header, translate).to_string();
    } else if revcomp == &1 {
        return format!("{} [revcomp]", header).to_string();
    } else if translate > &0 {
        return format!("{} [translate({})]", header, translate).to_string();
    }

    header.to_string()
}
pub fn read_fasta_string(contents: &String) -> HashMap<String, FastaResult> {
    // Initialize
    let mut hdr: Vec<&str> = Vec::new();
    let mut seq: String = String::from("");
    let mut result: HashMap<String, FastaResult> = HashMap::new();

    // GO through lines
    let mut x = 1;
    for line in contents.split("\n") {
        if line.starts_with(">") {
            if hdr.len() > 0 {
                let hdr_id = format!("{}{}", hdr[0].to_string(), x);
                result.insert(
                    hdr_id,
                    FastaResult {
                        coord_start: if hdr.len() > 2 {
                            hdr[1].parse::<u16>().unwrap()
                        } else {
                            0
                        },
                        coord_end: if hdr.len() > 2 {
                            hdr[2].parse::<u16>().unwrap()
                        } else {
                            0
                        },
                        sequence: seq,
                    },
                );
                seq = String::from("");
                x += 1;
            }

            hdr = line
                .trim_start_matches(">")
                .split(" ")
                .into_iter()
                .map(|c| c)
                .collect();
        } else {
            seq += line.trim_end();
        }
    }

    // Add remaining
    if hdr.len() > 0 {
        let hdr_id = format!("{}{}", hdr[0].to_string(), x);
        result.insert(
            hdr_id,
            FastaResult {
                coord_start: if hdr.len() > 2 {
                    hdr[1].parse::<u16>().unwrap()
                } else {
                    0
                },
                coord_end: if hdr.len() > 2 {
                    hdr[2].parse::<u16>().unwrap()
                } else {
                    0
                },
                sequence: seq,
            },
        );
    }

    result
}

pub fn get_overlap_percent(source: Range<u16>, dest: Range<u16>, is_rev: bool) -> Option<f32> {
    // Check for valid ranges
    if source.start > source.end || dest.start > dest.end {
        return None;
    }

    // Check for no overlap
    if dest.start > source.end || dest.end < source.start {
        return None;
    }

    // Get overlap start
    let start = if dest.start > source.start {
        dest.start
    } else {
        source.start
    };

    // Get end overlap
    let end = if dest.end < source.end {
        dest.end
    } else {
        source.end
    };

    if start > end {
        return None;
    }

    // Get lengths
    let length_overlap = end - start;
    let length_source = if is_rev {
        source.end - source.start
    } else {
        dest.end - dest.start
    };

    // Get percent overlap
    let percent: f32 = length_overlap as f32 / length_source as f32;
    Some(percent)
}
