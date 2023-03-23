use crate::models::{Hit, OrfTranscript};
use crate::algorithms::orf;
use log::{info, warn};
use biotools::CONFIG;


pub fn generate(hit: &Hit, orf: &OrfTranscript) -> Option<OrfTranscript> {

    // Check config
    if !CONFIG.switch.extend_orf {
        return None;
    }

    // Try to generate ofr from complete sequence
    let ext_orf = match orf::generate(&hit, true) {
        Some(r) => r,
        None => {
            warn!("Did not receive valid extended orf for hmmsearch id# {}, gene {}, using original ovr.", hit.hmmsearch_id, hit.gene_id);
            return None;
        }
    };

    // Ensure extended contains initial orf
    if ext_orf.cdna_start > orf.cdna_start_transcript || ext_orf.cdna_end < orf.cdna_end_transcript {
        warn!("Extended orf does not consume initial for hmm search id# {}, gene {}, (ext coords: {}-{}, initial coords: {}-{}), reverting to initial orv.", 
            hit.hmmsearch_id, hit.gene_id, ext_orf.cdna_start, ext_orf.cdna_end, orf.cdna_start_transcript, orf.cdna_end_transcript);
        return None;
    }

    // Info message
    info!("Found valid extended orf for hmm search id# {}, gene {}, (ext coords: {}-{}, initial coords: {}-{}), using instead of initial orv.", 
        hit.hmmsearch_id, hit.gene_id, ext_orf.cdna_start, ext_orf.cdna_end, orf.cdna_start_transcript, orf.cdna_end_transcript);

    // Check for any overlap
    if ext_orf.cdna_start > orf.cdna_end_transcript || ext_orf.cdna_end < orf.cdna_start_transcript {
        warn!("Extended orf does not contain any overlap for hmm search id# {}, gene {}, reverting to initial orf.", hit.hmmsearch_id, hit.gene_id);
        return None;
    } else if ext_orf.aa_start_transcript > hit.ali_end || ext_orf.aa_end_transcript < hit.ali_start {
        warn!("Extended orf does not contain any overlap for hmm search id# {}, gene {}, reverting to initial orf.", hit.hmmsearch_id, hit.gene_id);
        return None;
    }

    // Get start overlap
    let overlap_start = if ext_orf.aa_start_transcript > hit.ali_start {
        ext_orf.aa_start_transcript
    } else {
        hit.ali_start
    };

    // Get end overlap
    let overlap_end = if ext_orf.aa_end_transcript < hit.ali_end {
        ext_orf.aa_end_transcript
    } else {
        hit.ali_end
    };

    // Check overlap
    let overlap_percent: f32 = (overlap_end as f32 - overlap_start as f32) / (ext_orf.aa_end_transcript as f32 - ext_orf.aa_start_transcript as f32);
    if overlap_percent < CONFIG.search.min_overlap {
        warn!("Orf only overlaps extended orf by {} percent on hmm search id# {}, gene {}, reverting to initial orv", overlap_percent, hit.hmmsearch_id, hit.gene_id);
        return None;
    }

    Some(ext_orf)
}


