use std::collections::HashMap;
use std::str;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref NT_TABLE: HashMap<String, String> = instantiate_nt_table();
}

pub fn translate(sequence: &String) -> String {

    // Translate
    let chunks = sequence.as_bytes().chunks(3).map(|buf| unsafe { str::from_utf8_unchecked(buf) }).collect::<Vec<&str>>();

    let mut res: String = "".to_owned();
    for c in chunks {
        let s = match NT_TABLE.get(&c.to_string()) {
            Some(r) => r,
            None => continue
        };
        res.push_str(s.as_str());
    }




    res.to_string()
}

pub fn instantiate_nt_table() -> HashMap<String, String> {

    // Define nt table
    let mut nt_table: HashMap<String, String> = HashMap::new();
    nt_table.insert("TTT".to_string(), "F".to_string());
    nt_table.insert("TCT".to_string(), "S".to_string());
    nt_table.insert("TAT".to_string(), "Y".to_string());
    nt_table.insert("TGT".to_string(), "C".to_string());
    nt_table.insert("TTC".to_string(), "F".to_string());
    nt_table.insert("TCC".to_string(), "S".to_string());
    nt_table.insert("TAC".to_string(), "Y".to_string());
    nt_table.insert("TGC".to_string(), "C".to_string());
    nt_table.insert("TTA".to_string(), "L".to_string());
    nt_table.insert("TCA".to_string(), "S".to_string());
    nt_table.insert("TAA".to_string(), "*".to_string());
    nt_table.insert("TGA".to_string(), "*".to_string());
    nt_table.insert("TTG".to_string(), "L".to_string());
    nt_table.insert("TCG".to_string(), "S".to_string());
    nt_table.insert("TAG".to_string(), "*".to_string());
    nt_table.insert("TGG".to_string(), "W".to_string());
    nt_table.insert("CTT".to_string(), "L".to_string());
    nt_table.insert("CCT".to_string(), "P".to_string());
    nt_table.insert("CAT".to_string(), "H".to_string());
    nt_table.insert("CGT".to_string(), "R".to_string());
    nt_table.insert("CTC".to_string(), "L".to_string());
    nt_table.insert("CCC".to_string(), "P".to_string());
    nt_table.insert("CAC".to_string(), "H".to_string());
    nt_table.insert("CGC".to_string(), "R".to_string());
    nt_table.insert("CTA".to_string(), "L".to_string());
    nt_table.insert("CCA".to_string(), "P".to_string());
    nt_table.insert("CAA".to_string(), "Q".to_string());
    nt_table.insert("CGA".to_string(), "R".to_string());
    nt_table.insert("CTG".to_string(), "L".to_string());
    nt_table.insert("CCG".to_string(), "P".to_string());
    nt_table.insert("CAG".to_string(), "Q".to_string());
    nt_table.insert("CGG".to_string(), "R".to_string());
    nt_table.insert("ATT".to_string(), "I".to_string());
    nt_table.insert("ACT".to_string(), "T".to_string());
    nt_table.insert("AAT".to_string(), "N".to_string());
    nt_table.insert("AGT".to_string(), "S".to_string());
    nt_table.insert("ATC".to_string(), "I".to_string());
    nt_table.insert("ACC".to_string(), "T".to_string());
    nt_table.insert("AAC".to_string(), "N".to_string());
    nt_table.insert("AGC".to_string(), "S".to_string());
    nt_table.insert("ATA".to_string(), "I".to_string());
    nt_table.insert("ACA".to_string(), "T".to_string());
    nt_table.insert("AAA".to_string(), "K".to_string());
    nt_table.insert("AGA".to_string(), "R".to_string());
    nt_table.insert("ATG".to_string(), "M".to_string());
    nt_table.insert("ACG".to_string(), "T".to_string());
    nt_table.insert("AAG".to_string(), "K".to_string());
    nt_table.insert("AGG".to_string(), "R".to_string());
    nt_table.insert("GTT".to_string(), "V".to_string());
    nt_table.insert("GCT".to_string(), "A".to_string());
    nt_table.insert("GAT".to_string(), "D".to_string());
    nt_table.insert("GGT".to_string(), "G".to_string());
    nt_table.insert("GTC".to_string(), "V".to_string());
    nt_table.insert("GCC".to_string(), "A".to_string());
    nt_table.insert("GAC".to_string(), "D".to_string());
    nt_table.insert("GGC".to_string(), "G".to_string());
    nt_table.insert("GTA".to_string(), "V".to_string());
    nt_table.insert("GCA".to_string(), "A".to_string());
    nt_table.insert("GAA".to_string(), "E".to_string());
    nt_table.insert("GGA".to_string(), "G".to_string());
    nt_table.insert("GTG".to_string(), "V".to_string());
    nt_table.insert("GCG".to_string(), "A".to_string());
    nt_table.insert("GAG".to_string(), "E".to_string());
    nt_table.insert("GGG".to_string(), "G".to_string());

    // Return
    nt_table
}


