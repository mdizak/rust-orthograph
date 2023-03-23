use crate::algorithms::env_pseudo_master::EnvCandidate;
use crate::algorithms::frameshift_correction::OrfResult;
use crate::algorithms::hmm_overlap::HmmDiscard;
use crate::models::{Hit, HmmSearch};
use crate::reporter::ReporterKit;
use biotools::db::sqlite::TBL_HITS;
use biotools::CONFIG;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

pub struct Stats {
    nrh_by_gene: HashMap<String, u16>,
    discard_non_orf: u32,
    discard_hmm_overlap: u32,
    discard_env_pseudo_master: u32,
    discard_env_overlap: u32,
    pub discards: Vec<u32>,
    brh_fh: File,
    nolap_fh: File,
    sum_fh: File,
    filter_fh: File,
    report_fh: File,
}

impl Stats {
    pub fn new() -> Self {
        // Open files
        let brh_fh = biotools::io::open_file(format!(
            "{}/best-reciprocal-hits.txt",
            CONFIG.report.output_dir
        ));
        let nolap_fh = biotools::io::open_file(format!(
            "{}/non-overlapping-best-reciprocal-hits.txt",
            CONFIG.report.output_dir
        ));
        let sum_fh = biotools::io::open_file(format!("{}/summary.txt", CONFIG.report.output_dir));
        let report_fh = biotools::io::open_file(format!("{}/report.txt", CONFIG.report.output_dir));
        let filter_fh =
            biotools::io::open_file(format!("{}/filtered-hits.txt", CONFIG.report.output_dir));

        Self {
            nrh_by_gene: HashMap::new(),
            discard_non_orf: 0,
            discard_hmm_overlap: 0,
            discard_env_pseudo_master: 0,
            discard_env_overlap: 0,
            discards: Vec::new(),
            brh_fh: brh_fh,
            nolap_fh: nolap_fh,
            sum_fh: sum_fh,
            filter_fh: filter_fh,
            report_fh: report_fh,
        }
    }

    pub fn write_brh(&mut self, hit: &Hit) {
        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            hit.gene_id,
            hit.header_full,
            hit.ali_start,
            hit.ali_end,
            hit.score,
            hit.evalue,
            hit.hmm_start,
            hit.hmm_end
        );

        self.brh_fh
            .write_all(&line.as_bytes())
            .expect("Unable to write to best-recipocal-hits.txt file");
    }

    pub fn write_nolap(&mut self, hit: &Hit) {
        let line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            hit.gene_id,
            hit.header_full,
            hit.ali_start,
            hit.ali_end,
            hit.score,
            hit.evalue,
            hit.hmm_start,
            hit.hmm_end
        );

        self.nolap_fh
            .write_all(&line.as_bytes())
            .expect("Unable to write to non-overlapping-best-hits.txt file");
    }

    pub fn write_filtered_hit(
        &mut self,
        gene_id: &String,
        header: &String,
        revcomp: &u8,
        translate: &u8,
        reason: &str,
    ) {
        // Format line
        let tmp_header = biotools::format_header(&header, &revcomp, &translate);
        let line = format!("{},{},{}\n", gene_id, tmp_header, reason);

        self.filter_fh
            .write_all(&line.as_bytes())
            .expect("Unable to write to filtered-hits.txt file");
    }

    pub fn add_non_reciprocal_hit(&mut self, cand: &HmmSearch) {
        // Update counter
        let gene = format!("{}", cand.gene_id);
        *self.nrh_by_gene.entry(gene).or_insert(0) += 1;

        // Write to filtered hits file
        self.write_filtered_hit(&cand.gene_id, &cand.header, &0, &0, "non-reciprocal");
    }

    pub fn delete_hit(&mut self, kit: &ReporterKit, hit_id: &u32) {
        // Delete from in-memory db
        let sql = format!("DELETE FROM {} WHERE id = {}", *TBL_HITS, hit_id);
        match kit.memdb.execute(&sql, []) {
            Ok(_r) => {}
            Err(e) => panic!(
                "Unable to delete hit from temporary hits table, id# {}, error: {}",
                hit_id, e
            ),
        };
    }

    pub fn discard_non_orf(&mut self, kit: &ReporterKit, orf: &OrfResult) {
        // Unwrap variables
        let gene_id = orf.gene_id.as_ref().unwrap();
        let header_base = orf.header_base.as_ref().unwrap();
        let revcomp = orf.revcomp.as_ref().unwrap();
        let translate = orf.translate.as_ref().unwrap();

        // Discard hit
        self.delete_hit(&kit, &orf.hit_id);
        self.write_filtered_hit(&gene_id, &header_base, &revcomp, &translate, "no-orf-found");
        self.discard_non_orf += 1;
    }

    pub fn discard_hmm_overlap(&mut self, kit: &ReporterKit, hit: &HmmDiscard) {
        self.delete_hit(&kit, &hit.hit_id);
        self.write_filtered_hit(
            &hit.gene_id,
            &hit.header_base,
            &hit.revcomp,
            &hit.translate,
            "hmm-overlap",
        );
        self.discard_hmm_overlap += 1;
    }

    pub fn discard_env_pseudo_master(&mut self, kit: &ReporterKit, cand: &EnvCandidate) {
        self.delete_hit(&kit, &cand.id);
        self.write_filtered_hit(
            &cand.gene_id,
            &cand.header_base,
            &cand.hdr_revcomp,
            &cand.hdr_translate,
            "env-pseudo-master",
        );
        self.discard_env_pseudo_master += 1;
    }

    pub fn discard_env_overlap(&mut self, kit: &ReporterKit, cand: &EnvCandidate) {
        self.delete_hit(&kit, &cand.id);
        self.write_filtered_hit(
            &cand.gene_id,
            &cand.header_base,
            &cand.hdr_revcomp,
            &cand.hdr_translate,
            "env-overlap",
        );
        self.discard_env_overlap += 1;
    }

    pub fn write_report(&mut self) {
        self.report_fh
            .write_all("\n-- Report --\n\n".as_bytes())
            .expect("Unable to write to report.txt file");
        self.report_fh
            .write_all(
                format!(
                    "Skipped Env Pseudo Master: {}\n",
                    self.discard_env_pseudo_master
                )
                .as_bytes(),
            )
            .expect("Unable to write to report.txt file");
        self.report_fh
            .write_all(format!("Skipped Env Overlap: {}\n", self.discard_env_overlap).as_bytes())
            .expect("Unable to write to report.txt file");
        self.report_fh
            .write_all(format!("Skipped Hmm Overlap: {}\n", self.discard_hmm_overlap).as_bytes())
            .expect("Unable to write to report.txt file");
        self.report_fh
            .write_all(format!("Skipped No ORF: {}\n", self.discard_non_orf).as_bytes())
            .expect("Unable to write to report.txt file");
    }
}
