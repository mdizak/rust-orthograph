
pub struct HmmSearch {
    pub gene_id: String,
    pub aaseq_id: u32,
    pub ntseq_id: u32,
    pub taxid: u16,
    pub hmm_id: u32,
    pub score: f32,
    pub digest: String,
    pub evalue: String,
    pub hmm_start: u16,
    pub hmm_end: u16,
    pub ali_start: u16,
    pub ali_end: u16,
    pub env_start: u16,
    pub env_end: u16,
    pub header: String,
    pub non_orf_sequence: String
}

pub struct Hit {
    pub id: u32,
    pub is_overlap: bool,
    pub hmmsearch_id: u32,
    pub taxid: u16,
    pub aaseq_id: u32,
    pub ntseq_id: u32,
    pub blast_target: u32,
    pub gene_id: String,
    pub score: f32,
    pub digest: String,
    pub evalue: String,
    pub hmm_start: u16,
    pub hmm_end: u16,
    pub ali_start: u16,
    pub ali_end: u16,
    pub env_start: u16,
    pub env_end: u16,
    pub blast_start: u16,
    pub blast_end: u16,
    pub header_base: String,
    pub header_full: String,
    pub header_revcomp: bool,
    pub header_translate: u8,
    pub non_orf_sequence: String,
    pub est_sequence: String,
    pub hmm_sequence: String,
    pub aa_sequence: String
}

pub struct OrfTranscript {
    pub hit_id: u32,
    pub cdna_start: u16,
    pub cdna_end: u16,
    pub aa_start: u16,
    pub aa_end: u16,
    pub cdna_start_transcript: u16,
    pub cdna_end_transcript: u16,
    pub aa_start_transcript: u16,
    pub aa_end_transcript: u16,
    pub aa_start_hmm: u16,
    pub aa_end_hmm: u16,
    pub translated_seq: String,
    pub cdna_seq: String
}



