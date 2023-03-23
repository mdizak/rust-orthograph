use crate::models::HmmSearch;
use crate::reporter::ReporterKit;
use biotools::CONFIG;
use log::{info, warn};

pub fn check(kit: &ReporterKit, candidate: &HmmSearch) -> Option<(u32, u16, u16)> {

    // Get blast results
    info!("Getting blast results for '{}' (hmm search id# {}, alignment score {})", candidate.header, candidate.hmm_id, candidate.score); 
    let blasts = match kit.db.get_hmm_blast_results(candidate.hmm_id) {
        Ok(res) => res,
        Err(e) => panic!("Unable to obtain blast results for hmm search id# {}, error: {}", candidate.hmm_id, e)
    };

    // Check for zero blasts
    if blasts.len() == 0 {
        warn!("No blast results found for '{}' (gene '{}', hmm search id# {}), skipping.", candidate.header, candidate.gene_id, candidate.hmm_id);
        return None;
    }

    // Initialize
    let mut taxa_count = Vec::new();
    let mut mismatches: u16 = 0;

    // GO through blast results
    for num in 0..(blasts.len() - 1) {
        let blast = &blasts[num];

        // Get ref taxon name
        let ref_taxon: String = kit.db.get_ref_taxon_name(&blast.target).unwrap();

        // Check if hit occurs in hmm
        if kit.aaseq_by_gene[&candidate.gene_id].contains(&blast.target) {
            info!("    Reciprocal hit {} ({}) used in {}!", blast.target, ref_taxon, candidate.gene_id);

            // Check reference taxa
            if !&kit.reference_taxa.contains(&ref_taxon) {
                info!("'{}' not in reference taxon list, skipping", ref_taxon);
                continue;

            // Check if not in strict search, hence 1 hit ie enough
            } else if !CONFIG.switch.strict_search {
                return Some((blast.target, blast.res_start, blast.res_end));
            }

            // Add taxa to count, if not already threre
            if !taxa_count.contains(&ref_taxon) {
                taxa_count.push(ref_taxon);
            }

            // Check if we have all taxa = matches under strict search
            if taxa_count.len() >= kit.reference_taxa.len() {
                return Some((blast.target, blast.res_start, blast.res_end));
            }

        // Check one ahead for same score
        //} else if blasts.len() >= (num+1) && blasts[num+1].score == blast.score {
            //info!("    reciprocal hit {} ({}) (gene: {}) not used in this HMM, but next one has same score, skipping this hit", blast.target, ref_taxon, candidate.gene_id);
            //continue;

        // Mismatch
        } else {
            mismatches += 1;
            warn!("    reciprocal hit {} ({}) not used in this HMM (mismatch #{})", blast.target, ref_taxon, mismatches );

            // Check for too many mismatches
            if mismatches > CONFIG.search.max_mismatches {
                warn!("    Too many mismatches, we don't trust this one anymore.");
                return None;
            }
        }

    }

    // Not reciprocal
    None
}


