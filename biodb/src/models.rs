extern crate serde;

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug)]
pub struct HmmSearch {
    pub taxid: u16,
    pub gene: String,
    pub header: String,
    pub score: f32,
    pub evalue: String,
    pub log_evalue: f32,
    pub env_start: u16,
    pub env_end: u16,
    pub ali_start: u16,
    pub ali_end: u16,
    pub hmm_start: u16,
    pub hmm_end: u16,
    pub seq_type: u8,
    pub blast: Vec<Blast>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Blast {
    pub taxid: u16,
    pub target: u32,
    pub score: f32,
    pub evalue: String,
    pub log_evalue: f32,
    pub blast_start: u16,
    pub blast_end: u16,
    pub header: String,
    pub seq_type: u8,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Sequence {
    pub seq_type: u8,
    pub sequence: String,
}
