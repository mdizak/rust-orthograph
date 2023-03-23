use crate::models::Hit;
use crate::reporter::ReporterKit;
use biotools::db::sqlite::TBL_HITS;
use nclist::NClist;
use std::ops::Range;

pub fn check(kit: &ReporterKit, hit: &Hit, coords: &Vec<Range<u16>>) -> bool {
    // Get overlaps
    let nc = match NClist::from_vec(coords.to_vec()) {
        Ok(res) => res,
        Err(e) => panic!("Unable to check region mapped before overlap, error: {}", e),
    };

    // Check for overlaps
    let r = hit.ali_start..hit.ali_end;
    if nc.count_overlaps(&r) > 0 {
        return true;
    }

    // Update database
    let sql = format!("UPDATE {} SET is_overlap = 0 WHERE id = ?", *TBL_HITS);
    let mut stmt = match kit.memdb.prepare(&sql) {
        Ok(res) => res,
        Err(e) => panic!("Unable to prepare SQL statement to update no_overlap on temporary hits table, error: {}", e)
    };

    // Execute sql
    match stmt.query([&hit.id]) {
        Ok(r) => r,
        Err(e) => panic!("Unable to execute SQL statement to update no_overlap on temporary hits table, error: {}", e)
    };

    // No overlap found
    false
}
