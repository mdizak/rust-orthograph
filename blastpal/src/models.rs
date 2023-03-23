pub struct HmmSearch {
    pub id: u32,
    pub gene_id: String,
    pub target: String,
    pub ali_start: u16,
    pub ali_end: u16,
    pub sequence: String,
}

pub struct Blast {
    pub query: String,
    pub target: u16,
    pub score: f32,
    pub evalue: f32,
    pub log_evalue: f32,
    pub res_start: u16,
    pub res_end: u16,
    pub hmmsearch_id: u32,
}
