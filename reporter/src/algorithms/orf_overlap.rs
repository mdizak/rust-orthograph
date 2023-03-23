use crate::models::OrfTranscript;
use nclist::NClist;
use std::ops::Range;

pub fn check(orf: &OrfTranscript, coords: &Vec<Range<u16>>, is_extended: bool) -> bool {

    // Ensure we have coords to check
    if coords.len() == 0 {
        return false;
    }

    // Get coords
    let r = if is_extended { 
        orf.cdna_start..orf.cdna_end
    } else { 
        orf.cdna_start_transcript..orf.cdna_end_transcript
    };

println!("Overlap Chk: {}", orf.hit_id);
println!("Coord: {} / {}", r.start, r.end);
for c in coords { println!("Chk: {} / {}", c.start, c.end); }

    // Get overlaps
    let nc = match NClist::from_vec(coords.to_vec()) {
        Ok(res) => res,
        Err(e) => panic!("Unable to obtain overlap coords, error: {}", e)
    };

    // Check for overlaps
    if nc.count_overlaps(&r) > 0 {
        return true;
    }

    // No overlap found
    false
}

