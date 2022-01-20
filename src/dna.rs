use lib3dmol::structures::Structure;
use std::collections::HashMap;
use super::qt::Quaternion;
use super::constants::{INTERFACE_CUTOFF2, MEMBRANE_PENALTY_SCORE};
use super::scoring::{Score, satisfied_restraints, membrane_intersection};

macro_rules! hashmap {
    ($( $key: expr => $val: expr ),*) => {{
         let mut map = ::std::collections::HashMap::new();
         $( map.insert($key, $val); )*
         map
    }}
}

const EPSILON: f64 = 4.0;
const FACTOR: f64 = 332.0;
const MAX_ES_CUTOFF: f64 = 1.0;
const MIN_ES_CUTOFF: f64 = -1.0;
const VDW_CUTOFF: f64 = 1.0;
const ELEC_DIST_CUTOFF: f64 = 30.0;
const ELEC_DIST_CUTOFF2: f64 = ELEC_DIST_CUTOFF*ELEC_DIST_CUTOFF;
const VDW_DIST_CUTOFF: f64 = 10.0;
const VDW_DIST_CUTOFF2: f64 = VDW_DIST_CUTOFF*VDW_DIST_CUTOFF;
const ELEC_MAX_CUTOFF: f64 = MAX_ES_CUTOFF*EPSILON/FACTOR;
const ELEC_MIN_CUTOFF: f64 = MIN_ES_CUTOFF*EPSILON/FACTOR;

pub fn atoms_in_residues(residue_name: &str) -> &'static [&'static str] {
       match residue_name {
        "ALA" => { &["N", "CA", "C", "O", "CB"] }
        "CYS" => { &["N", "CA", "C", "O", "CB", "SG"] }
        "ASP" => { &["N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"] }
        "GLU" => { &["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"] }
        "PHE" => { &["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"] }
        "GLY" => { &["N", "CA", "C", "O"] }
        "HIS" => { &["N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"] }
        "ILE" => { &["N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"] }
        "LYS" => { &["N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"] }
        "LEU" => { &["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"] }
        "MET" => { &["N", "CA", "C", "O", "CB", "CG", "SD", "CE"] }
        "ASN" => { &["N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"] }
        "PRO" => { &["N", "CA", "C", "O", "CB", "CG", "CD"] }
        "GLN" => { &["N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"] }
        "ARG" => { &["N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"] }
        "SER" => { &["N", "CA", "C", "O", "CB", "OG"] }
        "THR" => { &["N", "CA", "C", "O", "CB", "OG1", "CG2"] }
        "VAL" => { &["N", "CA", "C", "O", "CB", "CG1", "CG2"] }
        "TRP" => { &["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE2", "NE1", "CE3", "CZ3", "CH2", "CZ2"] }
        "TYR" => { &["N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"] }
        "MMB" => { &["BJ"] }
        _ => { panic!("Residue name not supported in DNA scoring function") }
    }
}

lazy_static! {
    static ref VDW_CHARGES: HashMap<&'static str, f64> = hashmap![
        "IP" => 0.00277, "HS" => 0.0157, "HP" => 0.0157, "Na" => 0.00277, "N*" => 0.17, "Li" => 0.0183, "HO" => 0.0,
        "Rb" => 0.00017, "HC" => 0.0157, "HA" => 0.015, "O3" => 0.21, "CQ" => 0.086, "C*" => 0.086, "NA" => 0.17,
        "NB" => 0.17, "NC" => 0.17, "O2" => 0.21, "I" => 0.4, "Br" => 0.32, "H" => 0.0157, "HW" => 0.0, "C0" => 0.459789,
        "K" => 0.000328, "CK" => 0.086, "Cs" => 8.06e-05, "C" => 0.086, "Cl" => 0.1, "CN" => 0.086, "CM" => 0.086,
        "F" => 0.061, "CC" => 0.086, "CB" => 0.086, "CA" => 0.086, "Zn" => 0.0125, "O" => 0.21, "N" => 0.17,
        "P" => 0.2, "S" => 0.25, "CR" => 0.086, "N2" => 0.17, "N3" => 0.17, "CW" => 0.086, "CV" => 0.086, "CT" => 0.1094,
        "MG" => 0.8947, "OH" => 0.2104, "H2" => 0.0157, "H3" => 0.0157, "H1" => 0.0157, "H4" => 0.015, "H5" => 0.015,
        "SH" => 0.25, "OW" => 0.152, "OS" => 0.17];

    static ref VDW_RADII: HashMap<&'static str, f64> = hashmap![
        "IP" => 1.868, "HS" => 0.6, "HP" => 1.1, "Na" => 1.868, "N*" => 1.824, "Li" => 1.137, "HO" => 0.0001,
        "Rb" => 2.956, "HC" => 1.487, "HA" => 1.459, "O3" => 1.6612, "CQ" => 1.908, "C*" => 1.908,
        "NA" => 1.824, "NB" => 1.824, "NC" => 1.824, "O2" => 1.6612, "I" => 2.35, "Br" => 2.22,
        "H" => 0.6, "HW" => 0.0001, "C0" => 1.7131, "K" => 2.658, "CK" => 1.908, "Cs" => 3.395, "C" => 1.908,
        "Cl" => 2.47, "CN" => 1.908, "CM" => 1.908, "F" => 1.75, "CC" => 1.908, "CB" => 1.908, "CA" => 1.908,
        "Zn" => 1.1, "O" => 1.6612, "N" => 1.824, "P" => 2.1, "S" => 2.0, "CR" => 1.908, "N2" => 1.824,
        "N3" => 1.875, "CW" => 1.908, "CV" => 1.908, "CT" => 1.908, "MG" => 0.7926, "OH" => 1.721, "H2" => 1.287,
        "H3" => 1.187, "H1" => 1.387, "H4" => 1.409, "H5" => 1.359, "SH" => 2.0, "OW" => 1.7683, "OS" => 1.6837];

    static ref RES_TO_TRANSLATE: HashMap<&'static str, &'static str> = hashmap![
        "HIS" => "HID", "THY" => "DT", "ADE" => "DA", "CYT" => "DC", "GUA" => "DG"];

    static ref AMBER_TYPES: HashMap<&'static str, &'static str> = hashmap![
        "ALA-C" => "C", "ALA-CA" => "CT", "ALA-CB" => "CT", "ALA-H" => "H", "ALA-HA" => "H1", "ALA-HB1" => "HC", "ALA-HB2" => "HC", "ALA-HB3" => "HC", "ALA-N" => "N", "ALA-O" => "O",
        "ARG-C" => "C", "ARG-CA" => "CT", "ARG-CB" => "CT", "ARG-CD" => "CT", "ARG-CG" => "CT", "ARG-CZ" => "CA", "ARG-H" => "H", "ARG-HA" => "H1", "ARG-HB2" => "HC", "ARG-HB3" => "HC", "ARG-HD2" => "H1", "ARG-HD3" => "H1", "ARG-HE" => "H", "ARG-HG2" => "HC", "ARG-HG3" => "HC", "ARG-HH11" => "H", "ARG-HH12" => "H", "ARG-HH21" => "H", "ARG-HH22" => "H", "ARG-N" => "N", "ARG-NE" => "N2", "ARG-NH1" => "N2", "ARG-NH2" => "N2", "ARG-O" => "O",
        "ASH-C" => "C", "ASH-CA" => "CT", "ASH-CB" => "CT", "ASH-CG" => "C", "ASH-H" => "H", "ASH-HA" => "H1", "ASH-HB2" => "HC", "ASH-HB3" => "HC", "ASH-HD2" => "HO", "ASH-N" => "N", "ASH-O" => "O", "ASH-OD1" => "O", "ASH-OD2" => "OH",
        "ASN-C" => "C", "ASN-CA" => "CT", "ASN-CB" => "CT", "ASN-CG" => "C", "ASN-H" => "H", "ASN-HA" => "H1", "ASN-HB2" => "HC", "ASN-HB3" => "HC", "ASN-HD21" => "H", "ASN-HD22" => "H", "ASN-N" => "N", "ASN-ND2" => "N", "ASN-O" => "O", "ASN-OD1" => "O",
        "ASP-C" => "C", "ASP-CA" => "CT", "ASP-CB" => "CT", "ASP-CG" => "C", "ASP-H" => "H", "ASP-HA" => "H1", "ASP-HB2" => "HC", "ASP-HB3" => "HC", "ASP-N" => "N", "ASP-O" => "O", "ASP-OD1" => "O2", "ASP-OD2" => "O2",
        "CYM-C" => "C", "CYM-CA" => "CT", "CYM-CB" => "CT", "CYM-HA" => "H1", "CYM-HB2" => "H1", "CYM-HB3" => "H1", "CYM-HN" => "H", "CYM-N" => "N", "CYM-O" => "O", "CYM-SG" => "SH",
        "CYS-C" => "C", "CYS-CA" => "CT", "CYS-CB" => "CT", "CYS-H" => "H", "CYS-HA" => "H1", "CYS-HB2" => "H1", "CYS-HB3" => "H1", "CYS-HG" => "HS", "CYS-N" => "N", "CYS-O" => "O", "CYS-SG" => "SH",
        "CYX-C" => "C", "CYX-CA" => "CT", "CYX-CB" => "CT", "CYX-H" => "H", "CYX-HA" => "H1", "CYX-HB2" => "H1", "CYX-HB3" => "H1", "CYX-N" => "N", "CYX-O" => "O", "CYX-SG" => "S",
        "DA-C1'" => "CT", "DA-C2" => "CQ", "DA-C2'" => "CT", "DA-C3'" => "CT", "DA-C4" => "CB", "DA-C4'" => "CT", "DA-C5" => "CB", "DA-C5'" => "CT", "DA-C6" => "CA", "DA-C8" => "CK", "DA-H1'" => "H2", "DA-H2" => "H5", "DA-H2'1" => "HC", "DA-H2'2" => "HC", "DA-H3'" => "H1", "DA-H4'" => "H1", "DA-H5'1" => "H1", "DA-H5'2" => "H1", "DA-H61" => "H", "DA-H62" => "H", "DA-H8" => "H5", "DA-N1" => "NC", "DA-N3" => "NC", "DA-N6" => "N2", "DA-N7" => "NB", "DA-N9" => "N*", "DA-O1P" => "O2", "DA-O2P" => "O2", "DA-O3'" => "OS", "DA-O4'" => "OS", "DA-O5'" => "OS", "DA-P" => "P",
        "DA3-C1'" => "CT", "DA3-C2" => "CQ", "DA3-C2'" => "CT", "DA3-C3'" => "CT", "DA3-C4" => "CB", "DA3-C4'" => "CT", "DA3-C5" => "CB", "DA3-C5'" => "CT", "DA3-C6" => "CA", "DA3-C8" => "CK", "DA3-H1'" => "H2", "DA3-H2" => "H5", "DA3-H2'1" => "HC", "DA3-H2'2" => "HC", "DA3-H3'" => "H1", "DA3-H3T" => "HO", "DA3-H4'" => "H1", "DA3-H5'1" => "H1", "DA3-H5'2" => "H1", "DA3-H61" => "H", "DA3-H62" => "H", "DA3-H8" => "H5", "DA3-N1" => "NC", "DA3-N3" => "NC", "DA3-N6" => "N2", "DA3-N7" => "NB", "DA3-N9" => "N*", "DA3-O1P" => "O2", "DA3-O2P" => "O2", "DA3-O3'" => "OH", "DA3-O4'" => "OS", "DA3-O5'" => "OS", "DA3-P" => "P",
        "DA5-C1'" => "CT", "DA5-C2" => "CQ", "DA5-C2'" => "CT", "DA5-C3'" => "CT", "DA5-C4" => "CB", "DA5-C4'" => "CT", "DA5-C5" => "CB", "DA5-C5'" => "CT", "DA5-C6" => "CA", "DA5-C8" => "CK", "DA5-H1'" => "H2", "DA5-H2" => "H5", "DA5-H2'1" => "HC", "DA5-H2'2" => "HC", "DA5-H3'" => "H1", "DA5-H4'" => "H1", "DA5-H5'1" => "H1", "DA5-H5'2" => "H1", "DA5-H5T" => "HO", "DA5-H61" => "H", "DA5-H62" => "H", "DA5-H8" => "H5", "DA5-N1" => "NC", "DA5-N3" => "NC", "DA5-N6" => "N2", "DA5-N7" => "NB", "DA5-N9" => "N*", "DA5-O3'" => "OS", "DA5-O4'" => "OS", "DA5-O5'" => "OH",
        "DAN-C1'" => "CT", "DAN-C2" => "CQ", "DAN-C2'" => "CT", "DAN-C3'" => "CT", "DAN-C4" => "CB", "DAN-C4'" => "CT", "DAN-C5" => "CB", "DAN-C5'" => "CT", "DAN-C6" => "CA", "DAN-C8" => "CK", "DAN-H1'" => "H2", "DAN-H2" => "H5", "DAN-H2'1" => "HC", "DAN-H2'2" => "HC", "DAN-H3'" => "H1", "DAN-H3T" => "HO", "DAN-H4'" => "H1", "DAN-H5'1" => "H1", "DAN-H5'2" => "H1", "DAN-H5T" => "HO", "DAN-H61" => "H", "DAN-H62" => "H", "DAN-H8" => "H5", "DAN-N1" => "NC", "DAN-N3" => "NC", "DAN-N6" => "N2", "DAN-N7" => "NB", "DAN-N9" => "N*", "DAN-O3'" => "OH", "DAN-O4'" => "OS", "DAN-O5'" => "OH",
        "DC-C1'" => "CT", "DC-C2" => "C", "DC-C2'" => "CT", "DC-C3'" => "CT", "DC-C4" => "CA", "DC-C4'" => "CT", "DC-C5" => "CM", "DC-C5'" => "CT", "DC-C6" => "CM", "DC-H1'" => "H2", "DC-H2'1" => "HC", "DC-H2'2" => "HC", "DC-H3'" => "H1", "DC-H4'" => "H1", "DC-H41" => "H", "DC-H42" => "H", "DC-H5" => "HA", "DC-H5'1" => "H1", "DC-H5'2" => "H1", "DC-H6" => "H4", "DC-N1" => "N*", "DC-N3" => "NC", "DC-N4" => "N2", "DC-O1P" => "O2", "DC-O2" => "O", "DC-O2P" => "O2", "DC-O3'" => "OS", "DC-O4'" => "OS", "DC-O5'" => "OS", "DC-P" => "P",
        "DC3-C1'" => "CT", "DC3-C2" => "C", "DC3-C2'" => "CT", "DC3-C3'" => "CT", "DC3-C4" => "CA", "DC3-C4'" => "CT", "DC3-C5" => "CM", "DC3-C5'" => "CT", "DC3-C6" => "CM", "DC3-H1'" => "H2", "DC3-H2'1" => "HC", "DC3-H2'2" => "HC", "DC3-H3'" => "H1", "DC3-H3T" => "HO", "DC3-H4'" => "H1", "DC3-H41" => "H", "DC3-H42" => "H", "DC3-H5" => "HA", "DC3-H5'1" => "H1", "DC3-H5'2" => "H1", "DC3-H6" => "H4", "DC3-N1" => "N*", "DC3-N3" => "NC", "DC3-N4" => "N2", "DC3-O1P" => "O2", "DC3-O2" => "O", "DC3-O2P" => "O2", "DC3-O3'" => "OH", "DC3-O4'" => "OS", "DC3-O5'" => "OS", "DC3-P" => "P",
        "DC5-C1'" => "CT", "DC5-C2" => "C", "DC5-C2'" => "CT", "DC5-C3'" => "CT", "DC5-C4" => "CA", "DC5-C4'" => "CT", "DC5-C5" => "CM", "DC5-C5'" => "CT", "DC5-C6" => "CM", "DC5-H1'" => "H2", "DC5-H2'1" => "HC", "DC5-H2'2" => "HC", "DC5-H3'" => "H1", "DC5-H4'" => "H1", "DC5-H41" => "H", "DC5-H42" => "H", "DC5-H5" => "HA", "DC5-H5'1" => "H1", "DC5-H5'2" => "H1", "DC5-H5T" => "HO", "DC5-H6" => "H4", "DC5-N1" => "N*", "DC5-N3" => "NC", "DC5-N4" => "N2", "DC5-O2" => "O", "DC5-O3'" => "OS", "DC5-O4'" => "OS", "DC5-O5'" => "OH",
        "DCN-C1'" => "CT", "DCN-C2" => "C", "DCN-C2'" => "CT", "DCN-C3'" => "CT", "DCN-C4" => "CA", "DCN-C4'" => "CT", "DCN-C5" => "CM", "DCN-C5'" => "CT", "DCN-C6" => "CM", "DCN-H1'" => "H2", "DCN-H2'1" => "HC", "DCN-H2'2" => "HC", "DCN-H3'" => "H1", "DCN-H3T" => "HO", "DCN-H4'" => "H1", "DCN-H41" => "H", "DCN-H42" => "H", "DCN-H5" => "HA", "DCN-H5'1" => "H1", "DCN-H5'2" => "H1", "DCN-H5T" => "HO", "DCN-H6" => "H4", "DCN-N1" => "N*", "DCN-N3" => "NC", "DCN-N4" => "N2", "DCN-O2" => "O", "DCN-O3'" => "OH", "DCN-O4'" => "OS", "DCN-O5'" => "OH",
        "DG-C1'" => "CT", "DG-C2" => "CA", "DG-C2'" => "CT", "DG-C3'" => "CT", "DG-C4" => "CB", "DG-C4'" => "CT", "DG-C5" => "CB", "DG-C5'" => "CT", "DG-C6" => "C", "DG-C8" => "CK", "DG-H1" => "H", "DG-H1'" => "H2", "DG-H2'1" => "HC", "DG-H2'2" => "HC", "DG-H21" => "H", "DG-H22" => "H", "DG-H3'" => "H1", "DG-H4'" => "H1", "DG-H5'1" => "H1", "DG-H5'2" => "H1", "DG-H8" => "H5", "DG-N1" => "NA", "DG-N2" => "N2", "DG-N3" => "NC", "DG-N7" => "NB", "DG-N9" => "N*", "DG-O1P" => "O2", "DG-O2P" => "O2", "DG-O3'" => "OS", "DG-O4'" => "OS", "DG-O5'" => "OS", "DG-O6" => "O", "DG-P" => "P",
        "DG3-C1'" => "CT", "DG3-C2" => "CA", "DG3-C2'" => "CT", "DG3-C3'" => "CT", "DG3-C4" => "CB", "DG3-C4'" => "CT", "DG3-C5" => "CB", "DG3-C5'" => "CT", "DG3-C6" => "C", "DG3-C8" => "CK", "DG3-H1" => "H", "DG3-H1'" => "H2", "DG3-H2'1" => "HC", "DG3-H2'2" => "HC", "DG3-H21" => "H", "DG3-H22" => "H", "DG3-H3'" => "H1", "DG3-H3T" => "HO", "DG3-H4'" => "H1", "DG3-H5'1" => "H1", "DG3-H5'2" => "H1", "DG3-H8" => "H5", "DG3-N1" => "NA", "DG3-N2" => "N2", "DG3-N3" => "NC", "DG3-N7" => "NB", "DG3-N9" => "N*", "DG3-O1P" => "O2", "DG3-O2P" => "O2", "DG3-O3'" => "OH", "DG3-O4'" => "OS", "DG3-O5'" => "OS", "DG3-O6" => "O", "DG3-P" => "P",
        "DG5-C1'" => "CT", "DG5-C2" => "CA", "DG5-C2'" => "CT", "DG5-C3'" => "CT", "DG5-C4" => "CB", "DG5-C4'" => "CT", "DG5-C5" => "CB", "DG5-C5'" => "CT", "DG5-C6" => "C", "DG5-C8" => "CK", "DG5-H1" => "H", "DG5-H1'" => "H2", "DG5-H2'1" => "HC", "DG5-H2'2" => "HC", "DG5-H21" => "H", "DG5-H22" => "H", "DG5-H3'" => "H1", "DG5-H4'" => "H1", "DG5-H5'1" => "H1", "DG5-H5'2" => "H1", "DG5-H5T" => "HO", "DG5-H8" => "H5", "DG5-N1" => "NA", "DG5-N2" => "N2", "DG5-N3" => "NC", "DG5-N7" => "NB", "DG5-N9" => "N*", "DG5-O3'" => "OS", "DG5-O4'" => "OS", "DG5-O5'" => "OH", "DG5-O6" => "O",
        "DGN-C1'" => "CT", "DGN-C2" => "CA", "DGN-C2'" => "CT", "DGN-C3'" => "CT", "DGN-C4" => "CB", "DGN-C4'" => "CT", "DGN-C5" => "CB", "DGN-C5'" => "CT", "DGN-C6" => "C", "DGN-C8" => "CK", "DGN-H1" => "H", "DGN-H1'" => "H2", "DGN-H2'1" => "HC", "DGN-H2'2" => "HC", "DGN-H21" => "H", "DGN-H22" => "H", "DGN-H3'" => "H1", "DGN-H3T" => "HO", "DGN-H4'" => "H1", "DGN-H5'1" => "H1", "DGN-H5'2" => "H1", "DGN-H5T" => "HO", "DGN-H8" => "H5", "DGN-N1" => "NA", "DGN-N2" => "N2", "DGN-N3" => "NC", "DGN-N7" => "NB", "DGN-N9" => "N*", "DGN-O3'" => "OH", "DGN-O4'" => "OS", "DGN-O5'" => "OH", "DGN-O6" => "O",
        "DT-C1'" => "CT", "DT-C2" => "C", "DT-C2'" => "CT", "DT-C3'" => "CT", "DT-C4" => "C", "DT-C4'" => "CT", "DT-C5" => "CM", "DT-C5'" => "CT", "DT-C6" => "CM", "DT-C7" => "CT", "DT-H1'" => "H2", "DT-H2'1" => "HC", "DT-H2'2" => "HC", "DT-H3" => "H", "DT-H3'" => "H1", "DT-H4'" => "H1", "DT-H5'1" => "H1", "DT-H5'2" => "H1", "DT-H6" => "H4", "DT-H71" => "HC", "DT-H72" => "HC", "DT-H73" => "HC", "DT-N1" => "N*", "DT-N3" => "NA", "DT-O1P" => "O2", "DT-O2" => "O", "DT-O2P" => "O2", "DT-O3'" => "OS", "DT-O4" => "O", "DT-O4'" => "OS", "DT-O5'" => "OS", "DT-P" => "P",
        "DT3-C1'" => "CT", "DT3-C2" => "C", "DT3-C2'" => "CT", "DT3-C3'" => "CT", "DT3-C4" => "C", "DT3-C4'" => "CT", "DT3-C5" => "CM", "DT3-C5'" => "CT", "DT3-C6" => "CM", "DT3-C7" => "CT", "DT3-H1'" => "H2", "DT3-H2'1" => "HC", "DT3-H2'2" => "HC", "DT3-H3" => "H", "DT3-H3'" => "H1", "DT3-H3T" => "HO", "DT3-H4'" => "H1", "DT3-H5'1" => "H1", "DT3-H5'2" => "H1", "DT3-H6" => "H4", "DT3-H71" => "HC", "DT3-H72" => "HC", "DT3-H73" => "HC", "DT3-N1" => "N*", "DT3-N3" => "NA", "DT3-O1P" => "O2", "DT3-O2" => "O", "DT3-O2P" => "O2", "DT3-O3'" => "OH", "DT3-O4" => "O", "DT3-O4'" => "OS", "DT3-O5'" => "OS", "DT3-P" => "P",
        "DT5-C1'" => "CT", "DT5-C2" => "C", "DT5-C2'" => "CT", "DT5-C3'" => "CT", "DT5-C4" => "C", "DT5-C4'" => "CT", "DT5-C5" => "CM", "DT5-C5'" => "CT", "DT5-C6" => "CM", "DT5-C7" => "CT", "DT5-H1'" => "H2", "DT5-H2'1" => "HC", "DT5-H2'2" => "HC", "DT5-H3" => "H", "DT5-H3'" => "H1", "DT5-H4'" => "H1", "DT5-H5'1" => "H1", "DT5-H5'2" => "H1", "DT5-H5T" => "HO", "DT5-H6" => "H4", "DT5-H71" => "HC", "DT5-H72" => "HC", "DT5-H73" => "HC", "DT5-N1" => "N*", "DT5-N3" => "NA", "DT5-O2" => "O", "DT5-O3'" => "OS", "DT5-O4" => "O", "DT5-O4'" => "OS", "DT5-O5'" => "OH",
        "DTN-C1'" => "CT", "DTN-C2" => "C", "DTN-C2'" => "CT", "DTN-C3'" => "CT", "DTN-C4" => "C", "DTN-C4'" => "CT", "DTN-C5" => "CM", "DTN-C5'" => "CT", "DTN-C6" => "CM", "DTN-C7" => "CT", "DTN-H1'" => "H2", "DTN-H2'1" => "HC", "DTN-H2'2" => "HC", "DTN-H3" => "H", "DTN-H3'" => "H1", "DTN-H3T" => "HO", "DTN-H4'" => "H1", "DTN-H5'1" => "H1", "DTN-H5'2" => "H1", "DTN-H5T" => "HO", "DTN-H6" => "H4", "DTN-H71" => "HC", "DTN-H72" => "HC", "DTN-H73" => "HC", "DTN-N1" => "N*", "DTN-N3" => "NA", "DTN-O2" => "O", "DTN-O3'" => "OH", "DTN-O4" => "O", "DTN-O4'" => "OS", "DTN-O5'" => "OH",
        "GLH-C" => "C", "GLH-CA" => "CT", "GLH-CB" => "CT", "GLH-CD" => "C", "GLH-CG" => "CT", "GLH-H" => "H", "GLH-HA" => "H1", "GLH-HB2" => "HC", "GLH-HB3" => "HC", "GLH-HE2" => "HO", "GLH-HG2" => "HC", "GLH-HG3" => "HC", "GLH-N" => "N", "GLH-O" => "O", "GLH-OE1" => "O", "GLH-OE2" => "OH",
        "GLN-C" => "C", "GLN-CA" => "CT", "GLN-CB" => "CT", "GLN-CD" => "C", "GLN-CG" => "CT", "GLN-H" => "H", "GLN-HA" => "H1", "GLN-HB2" => "HC", "GLN-HB3" => "HC", "GLN-HE21" => "H", "GLN-HE22" => "H", "GLN-HG2" => "HC", "GLN-HG3" => "HC", "GLN-N" => "N", "GLN-NE2" => "N", "GLN-O" => "O", "GLN-OE1" => "O",
        "GLU-C" => "C", "GLU-CA" => "CT", "GLU-CB" => "CT", "GLU-CD" => "C", "GLU-CG" => "CT", "GLU-H" => "H", "GLU-HA" => "H1", "GLU-HB2" => "HC", "GLU-HB3" => "HC", "GLU-HG2" => "HC", "GLU-HG3" => "HC", "GLU-N" => "N", "GLU-O" => "O", "GLU-OE1" => "O2", "GLU-OE2" => "O2",
        "GLY-C" => "C", "GLY-CA" => "CT", "GLY-H" => "H", "GLY-HA2" => "H1", "GLY-HA3" => "H1", "GLY-N" => "N", "GLY-O" => "O",
        "HID-C" => "C", "HID-CA" => "CT", "HID-CB" => "CT", "HID-CD2" => "CV", "HID-CE1" => "CR", "HID-CG" => "CC", "HID-H" => "H", "HID-HA" => "H1", "HID-HB2" => "HC", "HID-HB3" => "HC", "HID-HD1" => "H", "HID-HD2" => "H4", "HID-HE1" => "H5", "HID-N" => "N", "HID-ND1" => "NA", "HID-NE2" => "NB", "HID-O" => "O",
        "HIS-C" => "C", "HIS-CA" => "CT", "HIS-CB" => "CT", "HIS-CD2" => "CV", "HIS-CE1" => "CR", "HIS-CG" => "CC", "HIS-H" => "H", "HIS-HA" => "H1", "HIS-HB2" => "HC", "HIS-HB3" => "HC", "HIS-HD1" => "H", "HIS-HD2" => "H4", "HIS-HE1" => "H5", "HIS-HE2" => "H", "HIS-N" => "N", "HIS-ND1" => "NA", "HIS-NE2" => "NB", "HIS-O" => "O",
        "HIE-C" => "C", "HIE-CA" => "CT", "HIE-CB" => "CT", "HIE-CD2" => "CW", "HIE-CE1" => "CR", "HIE-CG" => "CC", "HIE-H" => "H", "HIE-HA" => "H1", "HIE-HB2" => "HC", "HIE-HB3" => "HC", "HIE-HD2" => "H4", "HIE-HE1" => "H5", "HIE-HE2" => "H", "HIE-N" => "N", "HIE-ND1" => "NB", "HIE-NE2" => "NA", "HIE-O" => "O",
        "HIP-C" => "C", "HIP-CA" => "CT", "HIP-CB" => "CT", "HIP-CD2" => "CW", "HIP-CE1" => "CR", "HIP-CG" => "CC", "HIP-H" => "H", "HIP-HA" => "H1", "HIP-HB2" => "HC", "HIP-HB3" => "HC", "HIP-HD1" => "H", "HIP-HD2" => "H4", "HIP-HE1" => "H5", "HIP-HE2" => "H", "HIP-N" => "N", "HIP-ND1" => "NA", "HIP-NE2" => "NA", "HIP-O" => "O",
        "ILE-C" => "C", "ILE-CA" => "CT", "ILE-CB" => "CT", "ILE-CD1" => "CT", "ILE-CG1" => "CT", "ILE-CG2" => "CT", "ILE-H" => "H", "ILE-HA" => "H1", "ILE-HB" => "HC", "ILE-HD11" => "HC", "ILE-HD12" => "HC", "ILE-HD13" => "HC", "ILE-HG12" => "HC", "ILE-HG13" => "HC", "ILE-HG21" => "HC", "ILE-HG22" => "HC", "ILE-HG23" => "HC", "ILE-N" => "N", "ILE-O" => "O",
        "LEU-C" => "C", "LEU-CA" => "CT", "LEU-CB" => "CT", "LEU-CD1" => "CT", "LEU-CD2" => "CT", "LEU-CG" => "CT", "LEU-H" => "H", "LEU-HA" => "H1", "LEU-HB2" => "HC", "LEU-HB3" => "HC", "LEU-HD11" => "HC", "LEU-HD12" => "HC", "LEU-HD13" => "HC", "LEU-HD21" => "HC", "LEU-HD22" => "HC", "LEU-HD23" => "HC", "LEU-HG" => "HC", "LEU-N" => "N", "LEU-O" => "O",
        "LYN-C" => "C", "LYN-CA" => "CT", "LYN-CB" => "CT", "LYN-CD" => "CT", "LYN-CE" => "CT", "LYN-CG" => "CT", "LYN-H" => "H", "LYN-HA" => "H1", "LYN-HB2" => "HC", "LYN-HB3" => "HC", "LYN-HD2" => "HC", "LYN-HD3" => "HC", "LYN-HE2" => "HP", "LYN-HE3" => "HP", "LYN-HG2" => "HC", "LYN-HG3" => "HC", "LYN-HZ2" => "H", "LYN-HZ3" => "H", "LYN-N" => "N", "LYN-NZ" => "N3", "LYN-O" => "O",
        "LYS-C" => "C", "LYS-CA" => "CT", "LYS-CB" => "CT", "LYS-CD" => "CT", "LYS-CE" => "CT", "LYS-CG" => "CT", "LYS-H" => "H", "LYS-HA" => "H1", "LYS-HB2" => "HC", "LYS-HB3" => "HC", "LYS-HD2" => "HC", "LYS-HD3" => "HC", "LYS-HE2" => "HP", "LYS-HE3" => "HP", "LYS-HG2" => "HC", "LYS-HG3" => "HC", "LYS-HZ1" => "H", "LYS-HZ2" => "H", "LYS-HZ3" => "H", "LYS-N" => "N", "LYS-NZ" => "N3", "LYS-O" => "O",
        "MET-C" => "C", "MET-CA" => "CT", "MET-CB" => "CT", "MET-CE" => "CT", "MET-CG" => "CT", "MET-H" => "H", "MET-HA" => "H1", "MET-HB2" => "HC", "MET-HB3" => "HC", "MET-HE1" => "H1", "MET-HE2" => "H1", "MET-HE3" => "H1", "MET-HG2" => "H1", "MET-HG3" => "H1", "MET-N" => "N", "MET-O" => "O", "MET-SD" => "S",
        "PHE-C" => "C", "PHE-CA" => "CT", "PHE-CB" => "CT", "PHE-CD1" => "CA", "PHE-CD2" => "CA", "PHE-CE1" => "CA", "PHE-CE2" => "CA", "PHE-CG" => "CA", "PHE-CZ" => "CA", "PHE-H" => "H", "PHE-HA" => "H1", "PHE-HB2" => "HC", "PHE-HB3" => "HC", "PHE-HD1" => "HA", "PHE-HD2" => "HA", "PHE-HE1" => "HA", "PHE-HE2" => "HA", "PHE-HZ" => "HA", "PHE-N" => "N", "PHE-O" => "O",
        "PRO-C" => "C", "PRO-CA" => "CT", "PRO-CB" => "CT", "PRO-CD" => "CT", "PRO-CG" => "CT", "PRO-HA" => "H1", "PRO-HB2" => "HC", "PRO-HB3" => "HC", "PRO-HD2" => "H1", "PRO-HD3" => "H1", "PRO-HG2" => "HC", "PRO-HG3" => "HC", "PRO-N" => "N", "PRO-O" => "O",
        "RA-C1'" => "CT", "RA-C2" => "CQ", "RA-C2'" => "CT", "RA-C3'" => "CT", "RA-C4" => "CB", "RA-C4'" => "CT", "RA-C5" => "CB", "RA-C5'" => "CT", "RA-C6" => "CA", "RA-C8" => "CK", "RA-H1'" => "H2", "RA-H2" => "H5", "RA-H2'1" => "H1", "RA-H3'" => "H1", "RA-H4'" => "H1", "RA-H5'1" => "H1", "RA-H5'2" => "H1", "RA-H61" => "H", "RA-H62" => "H", "RA-H8" => "H5", "RA-HO'2" => "HO", "RA-N1" => "NC", "RA-N3" => "NC", "RA-N6" => "N2", "RA-N7" => "NB", "RA-N9" => "N*", "RA-O1P" => "O2", "RA-O2'" => "OH", "RA-O2P" => "O2", "RA-O3'" => "OS", "RA-O4'" => "OS", "RA-O5'" => "OS", "RA-P" => "P",
        "RA3-C1'" => "CT", "RA3-C2" => "CQ", "RA3-C2'" => "CT", "RA3-C3'" => "CT", "RA3-C4" => "CB", "RA3-C4'" => "CT", "RA3-C5" => "CB", "RA3-C5'" => "CT", "RA3-C6" => "CA", "RA3-C8" => "CK", "RA3-H1'" => "H2", "RA3-H2" => "H5", "RA3-H2'1" => "H1", "RA3-H3'" => "H1", "RA3-H3T" => "HO", "RA3-H4'" => "H1", "RA3-H5'1" => "H1", "RA3-H5'2" => "H1", "RA3-H61" => "H", "RA3-H62" => "H", "RA3-H8" => "H5", "RA3-HO'2" => "HO", "RA3-N1" => "NC", "RA3-N3" => "NC", "RA3-N6" => "N2", "RA3-N7" => "NB", "RA3-N9" => "N*", "RA3-O1P" => "O2", "RA3-O2'" => "OH", "RA3-O2P" => "O2", "RA3-O3'" => "OH", "RA3-O4'" => "OS", "RA3-O5'" => "OS", "RA3-P" => "P",
        "RA5-C1'" => "CT", "RA5-C2" => "CQ", "RA5-C2'" => "CT", "RA5-C3'" => "CT", "RA5-C4" => "CB", "RA5-C4'" => "CT", "RA5-C5" => "CB", "RA5-C5'" => "CT", "RA5-C6" => "CA", "RA5-C8" => "CK", "RA5-H1'" => "H2", "RA5-H2" => "H5", "RA5-H2'1" => "H1", "RA5-H3'" => "H1", "RA5-H4'" => "H1", "RA5-H5'1" => "H1", "RA5-H5'2" => "H1", "RA5-H5T" => "HO", "RA5-H61" => "H", "RA5-H62" => "H", "RA5-H8" => "H5", "RA5-HO'2" => "HO", "RA5-N1" => "NC", "RA5-N3" => "NC", "RA5-N6" => "N2", "RA5-N7" => "NB", "RA5-N9" => "N*", "RA5-O2'" => "OH", "RA5-O3'" => "OS", "RA5-O4'" => "OS", "RA5-O5'" => "OH",
        "RAN-C1'" => "CT", "RAN-C2" => "CQ", "RAN-C2'" => "CT", "RAN-C3'" => "CT", "RAN-C4" => "CB", "RAN-C4'" => "CT", "RAN-C5" => "CB", "RAN-C5'" => "CT", "RAN-C6" => "CA", "RAN-C8" => "CK", "RAN-H1'" => "H2", "RAN-H2" => "H5", "RAN-H2'1" => "H1", "RAN-H3'" => "H1", "RAN-H3T" => "HO", "RAN-H4'" => "H1", "RAN-H5'1" => "H1", "RAN-H5'2" => "H1", "RAN-H5T" => "HO", "RAN-H61" => "H", "RAN-H62" => "H", "RAN-H8" => "H5", "RAN-HO'2" => "HO", "RAN-N1" => "NC", "RAN-N3" => "NC", "RAN-N6" => "N2", "RAN-N7" => "NB", "RAN-N9" => "N*", "RAN-O2'" => "OH", "RAN-O3'" => "OH", "RAN-O4'" => "OS", "RAN-O5'" => "OH",
        "RC-C1'" => "CT", "RC-C2" => "C", "RC-C2'" => "CT", "RC-C3'" => "CT", "RC-C4" => "CA", "RC-C4'" => "CT", "RC-C5" => "CM", "RC-C5'" => "CT", "RC-C6" => "CM", "RC-H1'" => "H2", "RC-H2'1" => "H1", "RC-H3'" => "H1", "RC-H4'" => "H1", "RC-H41" => "H", "RC-H42" => "H", "RC-H5" => "HA", "RC-H5'1" => "H1", "RC-H5'2" => "H1", "RC-H6" => "H4", "RC-HO'2" => "HO", "RC-N1" => "N*", "RC-N3" => "NC", "RC-N4" => "N2", "RC-O1P" => "O2", "RC-O2" => "O", "RC-O2'" => "OH", "RC-O2P" => "O2", "RC-O3'" => "OS", "RC-O4'" => "OS", "RC-O5'" => "OS", "RC-P" => "P",
        "RC3-C1'" => "CT", "RC3-C2" => "C", "RC3-C2'" => "CT", "RC3-C3'" => "CT", "RC3-C4" => "CA", "RC3-C4'" => "CT", "RC3-C5" => "CM", "RC3-C5'" => "CT", "RC3-C6" => "CM", "RC3-H1'" => "H2", "RC3-H2'1" => "H1", "RC3-H3'" => "H1", "RC3-H3T" => "HO", "RC3-H4'" => "H1", "RC3-H41" => "H", "RC3-H42" => "H", "RC3-H5" => "HA", "RC3-H5'1" => "H1", "RC3-H5'2" => "H1", "RC3-H6" => "H4", "RC3-HO'2" => "HO", "RC3-N1" => "N*", "RC3-N3" => "NC", "RC3-N4" => "N2", "RC3-O1P" => "O2", "RC3-O2" => "O", "RC3-O2'" => "OH", "RC3-O2P" => "O2", "RC3-O3'" => "OH", "RC3-O4'" => "OS", "RC3-O5'" => "OS", "RC3-P" => "P",
        "RC5-C1'" => "CT", "RC5-C2" => "C", "RC5-C2'" => "CT", "RC5-C3'" => "CT", "RC5-C4" => "CA", "RC5-C4'" => "CT", "RC5-C5" => "CM", "RC5-C5'" => "CT", "RC5-C6" => "CM", "RC5-H1'" => "H2", "RC5-H2'1" => "H1", "RC5-H3'" => "H1", "RC5-H4'" => "H1", "RC5-H41" => "H", "RC5-H42" => "H", "RC5-H5" => "HA", "RC5-H5'1" => "H1", "RC5-H5'2" => "H1", "RC5-H5T" => "HO", "RC5-H6" => "H4", "RC5-HO'2" => "HO", "RC5-N1" => "N*", "RC5-N3" => "NC", "RC5-N4" => "N2", "RC5-O2" => "O", "RC5-O2'" => "OH", "RC5-O3'" => "OS", "RC5-O4'" => "OS", "RC5-O5'" => "OH",
        "RCN-C1'" => "CT", "RCN-C2" => "C", "RCN-C2'" => "CT", "RCN-C3'" => "CT", "RCN-C4" => "CA", "RCN-C4'" => "CT", "RCN-C5" => "CM", "RCN-C5'" => "CT", "RCN-C6" => "CM", "RCN-H1'" => "H2", "RCN-H2'1" => "H1", "RCN-H3'" => "H1", "RCN-H3T" => "HO", "RCN-H4'" => "H1", "RCN-H41" => "H", "RCN-H42" => "H", "RCN-H5" => "HA", "RCN-H5'1" => "H1", "RCN-H5'2" => "H1", "RCN-H5T" => "HO", "RCN-H6" => "H4", "RCN-HO'2" => "HO", "RCN-N1" => "N*", "RCN-N3" => "NC", "RCN-N4" => "N2", "RCN-O2" => "O", "RCN-O2'" => "OH", "RCN-O3'" => "OH", "RCN-O4'" => "OS", "RCN-O5'" => "OH",
        "RG-C1'" => "CT", "RG-C2" => "CA", "RG-C2'" => "CT", "RG-C3'" => "CT", "RG-C4" => "CB", "RG-C4'" => "CT", "RG-C5" => "CB", "RG-C5'" => "CT", "RG-C6" => "C", "RG-C8" => "CK", "RG-H1" => "H", "RG-H1'" => "H2", "RG-H2'1" => "H1", "RG-H21" => "H", "RG-H22" => "H", "RG-H3'" => "H1", "RG-H4'" => "H1", "RG-H5'1" => "H1", "RG-H5'2" => "H1", "RG-H8" => "H5", "RG-HO'2" => "HO", "RG-N1" => "NA", "RG-N2" => "N2", "RG-N3" => "NC", "RG-N7" => "NB", "RG-N9" => "N*", "RG-O1P" => "O2", "RG-O2'" => "OH", "RG-O2P" => "O2", "RG-O3'" => "OS", "RG-O4'" => "OS", "RG-O5'" => "OS", "RG-O6" => "O", "RG-P" => "P",
        "RG3-C1'" => "CT", "RG3-C2" => "CA", "RG3-C2'" => "CT", "RG3-C3'" => "CT", "RG3-C4" => "CB", "RG3-C4'" => "CT", "RG3-C5" => "CB", "RG3-C5'" => "CT", "RG3-C6" => "C", "RG3-C8" => "CK", "RG3-H1" => "H", "RG3-H1'" => "H2", "RG3-H2'1" => "H1", "RG3-H21" => "H", "RG3-H22" => "H", "RG3-H3'" => "H1", "RG3-H3T" => "HO", "RG3-H4'" => "H1", "RG3-H5'1" => "H1", "RG3-H5'2" => "H1", "RG3-H8" => "H5", "RG3-HO'2" => "HO", "RG3-N1" => "NA", "RG3-N2" => "N2", "RG3-N3" => "NC", "RG3-N7" => "NB", "RG3-N9" => "N*", "RG3-O1P" => "O2", "RG3-O2'" => "OH", "RG3-O2P" => "O2", "RG3-O3'" => "OH", "RG3-O4'" => "OS", "RG3-O5'" => "OS", "RG3-O6" => "O", "RG3-P" => "P",
        "RG5-C1'" => "CT", "RG5-C2" => "CA", "RG5-C2'" => "CT", "RG5-C3'" => "CT", "RG5-C4" => "CB", "RG5-C4'" => "CT", "RG5-C5" => "CB", "RG5-C5'" => "CT", "RG5-C6" => "C", "RG5-C8" => "CK", "RG5-H1" => "H", "RG5-H1'" => "H2", "RG5-H2'1" => "H1", "RG5-H21" => "H", "RG5-H22" => "H", "RG5-H3'" => "H1", "RG5-H4'" => "H1", "RG5-H5'1" => "H1", "RG5-H5'2" => "H1", "RG5-H5T" => "HO", "RG5-H8" => "H5", "RG5-HO'2" => "HO", "RG5-N1" => "NA", "RG5-N2" => "N2", "RG5-N3" => "NC", "RG5-N7" => "NB", "RG5-N9" => "N*", "RG5-O2'" => "OH", "RG5-O3'" => "OS", "RG5-O4'" => "OS", "RG5-O5'" => "OH", "RG5-O6" => "O",
        "RGN-C1'" => "CT", "RGN-C2" => "CA", "RGN-C2'" => "CT", "RGN-C3'" => "CT", "RGN-C4" => "CB", "RGN-C4'" => "CT", "RGN-C5" => "CB", "RGN-C5'" => "CT", "RGN-C6" => "C", "RGN-C8" => "CK", "RGN-H1" => "H", "RGN-H1'" => "H2", "RGN-H2'1" => "H1", "RGN-H21" => "H", "RGN-H22" => "H", "RGN-H3'" => "H1", "RGN-H3T" => "HO", "RGN-H4'" => "H1", "RGN-H5'1" => "H1", "RGN-H5'2" => "H1", "RGN-H5T" => "HO", "RGN-H8" => "H5", "RGN-HO'2" => "HO", "RGN-N1" => "NA", "RGN-N2" => "N2", "RGN-N3" => "NC", "RGN-N7" => "NB", "RGN-N9" => "N*", "RGN-O2'" => "OH", "RGN-O3'" => "OH", "RGN-O4'" => "OS", "RGN-O5'" => "OH", "RGN-O6" => "O",
        "RU-C1'" => "CT", "RU-C2" => "C", "RU-C2'" => "CT", "RU-C3'" => "CT", "RU-C4" => "C", "RU-C4'" => "CT", "RU-C5" => "CM", "RU-C5'" => "CT", "RU-C6" => "CM", "RU-H1'" => "H2", "RU-H2'1" => "H1", "RU-H3" => "H", "RU-H3'" => "H1", "RU-H4'" => "H1", "RU-H5" => "HA", "RU-H5'1" => "H1", "RU-H5'2" => "H1", "RU-H6" => "H4", "RU-HO'2" => "HO", "RU-N1" => "N*", "RU-N3" => "NA", "RU-O1P" => "O2", "RU-O2" => "O", "RU-O2'" => "OH", "RU-O2P" => "O2", "RU-O3'" => "OS", "RU-O4" => "O", "RU-O4'" => "OS", "RU-O5'" => "OS", "RU-P" => "P",
        "RU3-C1'" => "CT", "RU3-C2" => "C", "RU3-C2'" => "CT", "RU3-C3'" => "CT", "RU3-C4" => "C", "RU3-C4'" => "CT", "RU3-C5" => "CM", "RU3-C5'" => "CT", "RU3-C6" => "CM", "RU3-H1'" => "H2", "RU3-H2'1" => "H1", "RU3-H3" => "H", "RU3-H3'" => "H1", "RU3-H3T" => "HO", "RU3-H4'" => "H1", "RU3-H5" => "HA", "RU3-H5'1" => "H1", "RU3-H5'2" => "H1", "RU3-H6" => "H4", "RU3-HO'2" => "HO", "RU3-N1" => "N*", "RU3-N3" => "NA", "RU3-O1P" => "O2", "RU3-O2" => "O", "RU3-O2'" => "OH", "RU3-O2P" => "O2", "RU3-O3'" => "OH", "RU3-O4" => "O", "RU3-O4'" => "OS", "RU3-O5'" => "OS", "RU3-P" => "P",
        "RU5-C1'" => "CT", "RU5-C2" => "C", "RU5-C2'" => "CT", "RU5-C3'" => "CT", "RU5-C4" => "C", "RU5-C4'" => "CT", "RU5-C5" => "CM", "RU5-C5'" => "CT", "RU5-C6" => "CM", "RU5-H1'" => "H2", "RU5-H2'1" => "H1", "RU5-H3" => "H", "RU5-H3'" => "H1", "RU5-H4'" => "H1", "RU5-H5" => "HA", "RU5-H5'1" => "H1", "RU5-H5'2" => "H1", "RU5-H5T" => "HO", "RU5-H6" => "H4", "RU5-HO'2" => "HO", "RU5-N1" => "N*", "RU5-N3" => "NA", "RU5-O2" => "O", "RU5-O2'" => "OH", "RU5-O3'" => "OS", "RU5-O4" => "O", "RU5-O4'" => "OS", "RU5-O5'" => "OH",
        "RUN-C1'" => "CT", "RUN-C2" => "C", "RUN-C2'" => "CT", "RUN-C3'" => "CT", "RUN-C4" => "C", "RUN-C4'" => "CT", "RUN-C5" => "CM", "RUN-C5'" => "CT", "RUN-C6" => "CM", "RUN-H1'" => "H2", "RUN-H2'1" => "H1", "RUN-H3" => "H", "RUN-H3'" => "H1", "RUN-H3T" => "HO", "RUN-H4'" => "H1", "RUN-H5" => "HA", "RUN-H5'1" => "H1", "RUN-H5'2" => "H1", "RUN-H5T" => "HO", "RUN-H6" => "H4", "RUN-HO'2" => "HO", "RUN-N1" => "N*", "RUN-N3" => "NA", "RUN-O2" => "O", "RUN-O2'" => "OH", "RUN-O3'" => "OH", "RUN-O4" => "O", "RUN-O4'" => "OS", "RUN-O5'" => "OH",
        "SER-C" => "C", "SER-CA" => "CT", "SER-CB" => "CT", "SER-H" => "H", "SER-HA" => "H1", "SER-HB2" => "H1", "SER-HB3" => "H1", "SER-HG" => "HO", "SER-N" => "N", "SER-O" => "O", "SER-OG" => "OH",
        "THR-C" => "C", "THR-CA" => "CT", "THR-CB" => "CT", "THR-CG2" => "CT", "THR-H" => "H", "THR-HA" => "H1", "THR-HB" => "H1", "THR-HG1" => "HO", "THR-HG21" => "HC", "THR-HG22" => "HC", "THR-HG23" => "HC", "THR-N" => "N", "THR-O" => "O", "THR-OG1" => "OH",
        "TRP-C" => "C", "TRP-CA" => "CT", "TRP-CB" => "CT", "TRP-CD1" => "CW", "TRP-CD2" => "CB", "TRP-CE2" => "CN", "TRP-CE3" => "CA", "TRP-CG" => "C*", "TRP-CH2" => "CA", "TRP-CZ2" => "CA", "TRP-CZ3" => "CA", "TRP-H" => "H", "TRP-HA" => "H1", "TRP-HB2" => "HC", "TRP-HB3" => "HC", "TRP-HD1" => "H4", "TRP-HE1" => "H", "TRP-HE3" => "HA", "TRP-HH2" => "HA", "TRP-HZ2" => "HA", "TRP-HZ3" => "HA", "TRP-N" => "N", "TRP-NE1" => "NA", "TRP-O" => "O",
        "TYR-C" => "C", "TYR-CA" => "CT", "TYR-CB" => "CT", "TYR-CD1" => "CA", "TYR-CD2" => "CA", "TYR-CE1" => "CA", "TYR-CE2" => "CA", "TYR-CG" => "CA", "TYR-CZ" => "C", "TYR-H" => "H", "TYR-HA" => "H1", "TYR-HB2" => "HC", "TYR-HB3" => "HC", "TYR-HD1" => "HA", "TYR-HD2" => "HA", "TYR-HE1" => "HA", "TYR-HE2" => "HA", "TYR-HH" => "HO", "TYR-N" => "N", "TYR-O" => "O", "TYR-OH" => "OH",
        "VAL-C" => "C", "VAL-CA" => "CT", "VAL-CB" => "CT", "VAL-CG1" => "CT", "VAL-CG2" => "CT", "VAL-H" => "H", "VAL-HA" => "H1", "VAL-HB" => "HC", "VAL-HG11" => "HC", "VAL-HG12" => "HC", "VAL-HG13" => "HC", "VAL-HG21" => "HC", "VAL-HG22" => "HC", "VAL-HG23" => "HC", "VAL-N" => "N", "VAL-O" => "O"];

    static ref ELE_CHARGES: HashMap<&'static str, f64> = hashmap![
        "ALA-C" => 0.5973, "ALA-CA" => 0.0337, "ALA-CB" => -0.1825, "ALA-H" => 0.2719, "ALA-HA" => 0.0823, "ALA-HB1" => 0.0603, "ALA-HB2" => 0.0603, "ALA-HB3" => 0.0603, "ALA-N" => -0.4157, "ALA-O" => -0.5679,
        "ARG-C" => 0.7341, "ARG-CA" => -0.2637, "ARG-CB" => -0.0007, "ARG-CD" => 0.0486, "ARG-CG" => 0.039, "ARG-CZ" => 0.8076, "ARG-H" => 0.2747, "ARG-HA" => 0.156, "ARG-HB2" => 0.0327, "ARG-HB3" => 0.0327, "ARG-HD2" => 0.0687, "ARG-HD3" => 0.0687, "ARG-HE" => 0.3456, "ARG-HG2" => 0.0285, "ARG-HG3" => 0.0285, "ARG-HH11" => 0.4478, "ARG-HH12" => 0.4478, "ARG-HH21" => 0.4478, "ARG-HH22" => 0.4478, "ARG-N" => -0.3479, "ARG-NE" => -0.5295, "ARG-NH1" => -0.8627, "ARG-NH2" => -0.8627, "ARG-O" => -0.5894,
        "ASH-C" => 0.5973, "ASH-CA" => 0.0341, "ASH-CB" => -0.0316, "ASH-CG" => 0.6462, "ASH-H" => 0.2719, "ASH-HA" => 0.0864, "ASH-HB2" => 0.0488, "ASH-HB3" => 0.0488, "ASH-HD2" => 0.4747, "ASH-N" => -0.4157, "ASH-O" => -0.5679, "ASH-OD1" => -0.5554, "ASH-OD2" => -0.6376,
        "ASN-C" => 0.5973, "ASN-CA" => 0.0143, "ASN-CB" => -0.2041, "ASN-CG" => 0.713, "ASN-H" => 0.2719, "ASN-HA" => 0.1048, "ASN-HB2" => 0.0797, "ASN-HB3" => 0.0797, "ASN-HD21" => 0.4196, "ASN-HD22" => 0.4196, "ASN-N" => -0.4157, "ASN-ND2" => -0.9191, "ASN-O" => -0.5679, "ASN-OD1" => -0.5931,
        "ASP-C" => 0.5366, "ASP-CA" => 0.0381, "ASP-CB" => -0.0303, "ASP-CG" => 0.7994, "ASP-H" => 0.2936, "ASP-HA" => 0.088, "ASP-HB2" => -0.0122, "ASP-HB3" => -0.0122, "ASP-N" => -0.5163, "ASP-O" => -0.5819, "ASP-OD1" => -0.8014, "ASP-OD2" => -0.8014,
        "CYM-C" => 0.5973, "CYM-CA" => -0.0351, "CYM-CB" => -0.2413, "CYM-HA" => 0.0508, "CYM-HB2" => 0.1122, "CYM-HB3" => 0.1122, "CYM-HN" => 0.2719, "CYM-N" => -0.4157, "CYM-O" => -0.5679, "CYM-SG" => -0.8844,
        "CYS-C" => 0.5973, "CYS-CA" => 0.0213, "CYS-CB" => -0.1231, "CYS-H" => 0.2719, "CYS-HA" => 0.1124, "CYS-HB2" => 0.1112, "CYS-HB3" => 0.1112, "CYS-HG" => 0.1933, "CYS-N" => -0.4157, "CYS-O" => -0.5679, "CYS-SG" => -0.3119,
        "CYX-C" => 0.5973, "CYX-CA" => 0.0429, "CYX-CB" => -0.079, "CYX-H" => 0.2719, "CYX-HA" => 0.0766, "CYX-HB2" => 0.091, "CYX-HB3" => 0.091, "CYX-N" => -0.4157, "CYX-O" => -0.5679, "CYX-SG" => -0.1081,
        "DA-C1'" => 0.0431, "DA-C2" => 0.5716, "DA-C2'" => -0.0854, "DA-C3'" => 0.0713, "DA-C4" => 0.38, "DA-C4'" => 0.1629, "DA-C5" => 0.0725, "DA-C5'" => -0.0069, "DA-C6" => 0.6897, "DA-C8" => 0.1607, "DA-H1'" => 0.1838, "DA-H2" => 0.0598, "DA-H2'1" => 0.0718, "DA-H2'2" => 0.0718, "DA-H3'" => 0.0985, "DA-H4'" => 0.1176, "DA-H5'1" => 0.0754, "DA-H5'2" => 0.0754, "DA-H61" => 0.4167, "DA-H62" => 0.4167, "DA-H8" => 0.1877, "DA-N1" => -0.7624, "DA-N3" => -0.7417, "DA-N6" => -0.9123, "DA-N7" => -0.6175, "DA-N9" => -0.0268, "DA-O1P" => -0.7761, "DA-O2P" => -0.7761, "DA-O3'" => -0.5232, "DA-O4'" => -0.3691, "DA-O5'" => -0.4954, "DA-P" => 1.1659,
        "DA3-C1'" => 0.0431, "DA3-C2" => 0.5716, "DA3-C2'" => -0.0854, "DA3-C3'" => 0.0713, "DA3-C4" => 0.38, "DA3-C4'" => 0.1629, "DA3-C5" => 0.0725, "DA3-C5'" => -0.0069, "DA3-C6" => 0.6897, "DA3-C8" => 0.1607, "DA3-H1'" => 0.1838, "DA3-H2" => 0.0598, "DA3-H2'1" => 0.0718, "DA3-H2'2" => 0.0718, "DA3-H3'" => 0.0985, "DA3-H3T" => 0.4396, "DA3-H4'" => 0.1176, "DA3-H5'1" => 0.0754, "DA3-H5'2" => 0.0754, "DA3-H61" => 0.4167, "DA3-H62" => 0.4167, "DA3-H8" => 0.1877, "DA3-N1" => -0.7624, "DA3-N3" => -0.7417, "DA3-N6" => -0.9123, "DA3-N7" => -0.6175, "DA3-N9" => -0.0268, "DA3-O1P" => -0.7761, "DA3-O2P" => -0.7761, "DA3-O3'" => -0.6549, "DA3-O4'" => -0.3691, "DA3-O5'" => -0.4954, "DA3-P" => 1.1659,
        "DA5-C1'" => 0.0431, "DA5-C2" => 0.5716, "DA5-C2'" => -0.0854, "DA5-C3'" => 0.0713, "DA5-C4" => 0.38, "DA5-C4'" => 0.1629, "DA5-C5" => 0.0725, "DA5-C5'" => -0.0069, "DA5-C6" => 0.6897, "DA5-C8" => 0.1607, "DA5-H1'" => 0.1838, "DA5-H2" => 0.0598, "DA5-H2'1" => 0.0718, "DA5-H2'2" => 0.0718, "DA5-H3'" => 0.0985, "DA5-H4'" => 0.1176, "DA5-H5'1" => 0.0754, "DA5-H5'2" => 0.0754, "DA5-H5T" => 0.4422, "DA5-H61" => 0.4167, "DA5-H62" => 0.4167, "DA5-H8" => 0.1877, "DA5-N1" => -0.7624, "DA5-N3" => -0.7417, "DA5-N6" => -0.9123, "DA5-N7" => -0.6175, "DA5-N9" => -0.0268, "DA5-O3'" => -0.5232, "DA5-O4'" => -0.3691, "DA5-O5'" => -0.6318,
        "DAN-C1'" => 0.0431, "DAN-C2" => 0.5716, "DAN-C2'" => -0.0854, "DAN-C3'" => 0.0713, "DAN-C4" => 0.38, "DAN-C4'" => 0.1629, "DAN-C5" => 0.0725, "DAN-C5'" => -0.0069, "DAN-C6" => 0.6897, "DAN-C8" => 0.1607, "DAN-H1'" => 0.1838, "DAN-H2" => 0.0598, "DAN-H2'1" => 0.0718, "DAN-H2'2" => 0.0718, "DAN-H3'" => 0.0985, "DAN-H3T" => 0.4396, "DAN-H4'" => 0.1176, "DAN-H5'1" => 0.0754, "DAN-H5'2" => 0.0754, "DAN-H5T" => 0.4422, "DAN-H61" => 0.4167, "DAN-H62" => 0.4167, "DAN-H8" => 0.1877, "DAN-N1" => -0.7624, "DAN-N3" => -0.7417, "DAN-N6" => -0.9123, "DAN-N7" => -0.6175, "DAN-N9" => -0.0268, "DAN-O3'" => -0.6549, "DAN-O4'" => -0.3691, "DAN-O5'" => -0.6318,
        "DC-C1'" => -0.0116, "DC-C2" => 0.7959, "DC-C2'" => -0.0854, "DC-C3'" => 0.0713, "DC-C4" => 0.8439, "DC-C4'" => 0.1629, "DC-C5" => -0.5222, "DC-C5'" => -0.0069, "DC-C6" => -0.0183, "DC-H1'" => 0.1963, "DC-H2'1" => 0.0718, "DC-H2'2" => 0.0718, "DC-H3'" => 0.0985, "DC-H4'" => 0.1176, "DC-H41" => 0.4314, "DC-H42" => 0.4314, "DC-H5" => 0.1863, "DC-H5'1" => 0.0754, "DC-H5'2" => 0.0754, "DC-H6" => 0.2293, "DC-N1" => -0.0339, "DC-N3" => -0.7748, "DC-N4" => -0.9773, "DC-O1P" => -0.7761, "DC-O2" => -0.6548, "DC-O2P" => -0.7761, "DC-O3'" => -0.5232, "DC-O4'" => -0.3691, "DC-O5'" => -0.4954, "DC-P" => 1.1659,
        "DC3-C1'" => -0.0116, "DC3-C2" => 0.7959, "DC3-C2'" => -0.0854, "DC3-C3'" => 0.0713, "DC3-C4" => 0.8439, "DC3-C4'" => 0.1629, "DC3-C5" => -0.5222, "DC3-C5'" => -0.0069, "DC3-C6" => -0.0183, "DC3-H1'" => 0.1963, "DC3-H2'1" => 0.0718, "DC3-H2'2" => 0.0718, "DC3-H3'" => 0.0985, "DC3-H3T" => 0.4396, "DC3-H4'" => 0.1176, "DC3-H41" => 0.4314, "DC3-H42" => 0.4314, "DC3-H5" => 0.1863, "DC3-H5'1" => 0.0754, "DC3-H5'2" => 0.0754, "DC3-H6" => 0.2293, "DC3-N1" => -0.0339, "DC3-N3" => -0.7748, "DC3-N4" => -0.9773, "DC3-O1P" => -0.7761, "DC3-O2" => -0.6548, "DC3-O2P" => -0.7761, "DC3-O3'" => -0.6549, "DC3-O4'" => -0.3691, "DC3-O5'" => -0.4954, "DC3-P" => 1.1659,
        "DC5-C1'" => -0.0116, "DC5-C2" => 0.7959, "DC5-C2'" => -0.0854, "DC5-C3'" => 0.0713, "DC5-C4" => 0.8439, "DC5-C4'" => 0.1629, "DC5-C5" => -0.5222, "DC5-C5'" => -0.0069, "DC5-C6" => -0.0183, "DC5-H1'" => 0.1963, "DC5-H2'1" => 0.0718, "DC5-H2'2" => 0.0718, "DC5-H3'" => 0.0985, "DC5-H4'" => 0.1176, "DC5-H41" => 0.4314, "DC5-H42" => 0.4314, "DC5-H5" => 0.1863, "DC5-H5'1" => 0.0754, "DC5-H5'2" => 0.0754, "DC5-H5T" => 0.4422, "DC5-H6" => 0.2293, "DC5-N1" => -0.0339, "DC5-N3" => -0.7748, "DC5-N4" => -0.9773, "DC5-O2" => -0.6548, "DC5-O3'" => -0.5232, "DC5-O4'" => -0.3691, "DC5-O5'" => -0.6318,
        "DCN-C1'" => -0.0116, "DCN-C2" => 0.7959, "DCN-C2'" => -0.0854, "DCN-C3'" => 0.0713, "DCN-C4" => 0.8439, "DCN-C4'" => 0.1629, "DCN-C5" => -0.5222, "DCN-C5'" => -0.0069, "DCN-C6" => -0.0183, "DCN-H1'" => 0.1963, "DCN-H2'1" => 0.0718, "DCN-H2'2" => 0.0718, "DCN-H3'" => 0.0985, "DCN-H3T" => 0.4396, "DCN-H4'" => 0.1176, "DCN-H41" => 0.4314, "DCN-H42" => 0.4314, "DCN-H5" => 0.1863, "DCN-H5'1" => 0.0754, "DCN-H5'2" => 0.0754, "DCN-H5T" => 0.4422, "DCN-H6" => 0.2293, "DCN-N1" => -0.0339, "DCN-N3" => -0.7748, "DCN-N4" => -0.9773, "DCN-O2" => -0.6548, "DCN-O3'" => -0.6549, "DCN-O4'" => -0.3691, "DCN-O5'" => -0.6318,
        "DG-C1'" => 0.0358, "DG-C2" => 0.7432, "DG-C2'" => -0.0854, "DG-C3'" => 0.0713, "DG-C4" => 0.1814, "DG-C4'" => 0.1629, "DG-C5" => 0.1991, "DG-C5'" => -0.0069, "DG-C6" => 0.4918, "DG-C8" => 0.0736, "DG-H1" => 0.352, "DG-H1'" => 0.1746, "DG-H2'1" => 0.0718, "DG-H2'2" => 0.0718, "DG-H21" => 0.4235, "DG-H22" => 0.4235, "DG-H3'" => 0.0985, "DG-H4'" => 0.1176, "DG-H5'1" => 0.0754, "DG-H5'2" => 0.0754, "DG-H8" => 0.1997, "DG-N1" => -0.5053, "DG-N2" => -0.923, "DG-N3" => -0.6636, "DG-N7" => -0.5725, "DG-N9" => 0.0577, "DG-O1P" => -0.7761, "DG-O2P" => -0.7761, "DG-O3'" => -0.5232, "DG-O4'" => -0.3691, "DG-O5'" => -0.4954, "DG-O6" => -0.5699, "DG-P" => 1.1659,
        "DG3-C1'" => 0.0358, "DG3-C2" => 0.7432, "DG3-C2'" => -0.0854, "DG3-C3'" => 0.0713, "DG3-C4" => 0.1814, "DG3-C4'" => 0.1629, "DG3-C5" => 0.1991, "DG3-C5'" => -0.0069, "DG3-C6" => 0.4918, "DG3-C8" => 0.0736, "DG3-H1" => 0.352, "DG3-H1'" => 0.1746, "DG3-H2'1" => 0.0718, "DG3-H2'2" => 0.0718, "DG3-H21" => 0.4235, "DG3-H22" => 0.4235, "DG3-H3'" => 0.0985, "DG3-H3T" => 0.4396, "DG3-H4'" => 0.1176, "DG3-H5'1" => 0.0754, "DG3-H5'2" => 0.0754, "DG3-H8" => 0.1997, "DG3-N1" => -0.5053, "DG3-N2" => -0.923, "DG3-N3" => -0.6636, "DG3-N7" => -0.5725, "DG3-N9" => 0.0577, "DG3-O1P" => -0.7761, "DG3-O2P" => -0.7761, "DG3-O3'" => -0.6549, "DG3-O4'" => -0.3691, "DG3-O5'" => -0.4954, "DG3-O6" => -0.5699, "DG3-P" => 1.1659,
        "DG5-C1'" => 0.0358, "DG5-C2" => 0.7432, "DG5-C2'" => -0.0854, "DG5-C3'" => 0.0713, "DG5-C4" => 0.1814, "DG5-C4'" => 0.1629, "DG5-C5" => 0.1991, "DG5-C5'" => -0.0069, "DG5-C6" => 0.4918, "DG5-C8" => 0.0736, "DG5-H1" => 0.352, "DG5-H1'" => 0.1746, "DG5-H2'1" => 0.0718, "DG5-H2'2" => 0.0718, "DG5-H21" => 0.4235, "DG5-H22" => 0.4235, "DG5-H3'" => 0.0985, "DG5-H4'" => 0.1176, "DG5-H5'1" => 0.0754, "DG5-H5'2" => 0.0754, "DG5-H5T" => 0.4422, "DG5-H8" => 0.1997, "DG5-N1" => -0.5053, "DG5-N2" => -0.923, "DG5-N3" => -0.6636, "DG5-N7" => -0.5725, "DG5-N9" => 0.0577, "DG5-O3'" => -0.5232, "DG5-O4'" => -0.3691, "DG5-O5'" => -0.6318, "DG5-O6" => -0.5699,
        "DGN-C1'" => 0.0358, "DGN-C2" => 0.7432, "DGN-C2'" => -0.0854, "DGN-C3'" => 0.0713, "DGN-C4" => 0.1814, "DGN-C4'" => 0.1629, "DGN-C5" => 0.1991, "DGN-C5'" => -0.0069, "DGN-C6" => 0.4918, "DGN-C8" => 0.0736, "DGN-H1" => 0.352, "DGN-H1'" => 0.1746, "DGN-H2'1" => 0.0718, "DGN-H2'2" => 0.0718, "DGN-H21" => 0.4235, "DGN-H22" => 0.4235, "DGN-H3'" => 0.0985, "DGN-H3T" => 0.4396, "DGN-H4'" => 0.1176, "DGN-H5'1" => 0.0754, "DGN-H5'2" => 0.0754, "DGN-H5T" => 0.4422, "DGN-H8" => 0.1997, "DGN-N1" => -0.5053, "DGN-N2" => -0.923, "DGN-N3" => -0.6636, "DGN-N7" => -0.5725, "DGN-N9" => 0.0577, "DGN-O3'" => -0.6549, "DGN-O4'" => -0.3691, "DGN-O5'" => -0.6318, "DGN-O6" => -0.5699,
        "DT-C1'" => 0.068, "DT-C2" => 0.5677, "DT-C2'" => -0.0854, "DT-C3'" => 0.0713, "DT-C4" => 0.5194, "DT-C4'" => 0.1629, "DT-C5" => 0.0025, "DT-C5'" => -0.0069, "DT-C6" => -0.2209, "DT-C7" => -0.2269, "DT-H1'" => 0.1804, "DT-H2'1" => 0.0718, "DT-H2'2" => 0.0718, "DT-H3" => 0.342, "DT-H3'" => 0.0985, "DT-H4'" => 0.1176, "DT-H5'1" => 0.0754, "DT-H5'2" => 0.0754, "DT-H6" => 0.2607, "DT-H71" => 0.077, "DT-H72" => 0.077, "DT-H73" => 0.077, "DT-N1" => -0.0239, "DT-N3" => -0.434, "DT-O1P" => -0.7761, "DT-O2" => -0.5881, "DT-O2P" => -0.7761, "DT-O3'" => -0.5232, "DT-O4" => -0.5563, "DT-O4'" => -0.3691, "DT-O5'" => -0.4954, "DT-P" => 1.1659,
        "DT3-C1'" => 0.068, "DT3-C2" => 0.5677, "DT3-C2'" => -0.0854, "DT3-C3'" => 0.0713, "DT3-C4" => 0.5194, "DT3-C4'" => 0.1629, "DT3-C5" => 0.0025, "DT3-C5'" => -0.0069, "DT3-C6" => -0.2209, "DT3-C7" => -0.2269, "DT3-H1'" => 0.1804, "DT3-H2'1" => 0.0718, "DT3-H2'2" => 0.0718, "DT3-H3" => 0.342, "DT3-H3'" => 0.0985, "DT3-H3T" => 0.4396, "DT3-H4'" => 0.1176, "DT3-H5'1" => 0.0754, "DT3-H5'2" => 0.0754, "DT3-H6" => 0.2607, "DT3-H71" => 0.077, "DT3-H72" => 0.077, "DT3-H73" => 0.077, "DT3-N1" => -0.0239, "DT3-N3" => -0.434, "DT3-O1P" => -0.7761, "DT3-O2" => -0.5881, "DT3-O2P" => -0.7761, "DT3-O3'" => -0.6549, "DT3-O4" => -0.5563, "DT3-O4'" => -0.3691, "DT3-O5'" => -0.4954, "DT3-P" => 1.1659,
        "DT5-C1'" => 0.068, "DT5-C2" => 0.5677, "DT5-C2'" => -0.0854, "DT5-C3'" => 0.0713, "DT5-C4" => 0.5194, "DT5-C4'" => 0.1629, "DT5-C5" => 0.0025, "DT5-C5'" => -0.0069, "DT5-C6" => -0.2209, "DT5-C7" => -0.2269, "DT5-H1'" => 0.1804, "DT5-H2'1" => 0.0718, "DT5-H2'2" => 0.0718, "DT5-H3" => 0.342, "DT5-H3'" => 0.0985, "DT5-H4'" => 0.1176, "DT5-H5'1" => 0.0754, "DT5-H5'2" => 0.0754, "DT5-H5T" => 0.4422, "DT5-H6" => 0.2607, "DT5-H71" => 0.077, "DT5-H72" => 0.077, "DT5-H73" => 0.077, "DT5-N1" => -0.0239, "DT5-N3" => -0.434, "DT5-O2" => -0.5881, "DT5-O3'" => -0.5232, "DT5-O4" => -0.5563, "DT5-O4'" => -0.3691, "DT5-O5'" => -0.6318,
        "DTN-C1'" => 0.068, "DTN-C2" => 0.5677, "DTN-C2'" => -0.0854, "DTN-C3'" => 0.0713, "DTN-C4" => 0.5194, "DTN-C4'" => 0.1629, "DTN-C5" => 0.0025, "DTN-C5'" => -0.0069, "DTN-C6" => -0.2209, "DTN-C7" => -0.2269, "DTN-H1'" => 0.1804, "DTN-H2'1" => 0.0718, "DTN-H2'2" => 0.0718, "DTN-H3" => 0.342, "DTN-H3'" => 0.0985, "DTN-H3T" => 0.4396, "DTN-H4'" => 0.1176, "DTN-H5'1" => 0.0754, "DTN-H5'2" => 0.0754, "DTN-H5T" => 0.4422, "DTN-H6" => 0.2607, "DTN-H71" => 0.077, "DTN-H72" => 0.077, "DTN-H73" => 0.077, "DTN-N1" => -0.0239, "DTN-N3" => -0.434, "DTN-O2" => -0.5881, "DTN-O3'" => -0.6549, "DTN-O4" => -0.5563, "DTN-O4'" => -0.3691, "DTN-O5'" => -0.6318,
        "GLH-C" => 0.5973, "GLH-CA" => 0.0145, "GLH-CB" => -0.0071, "GLH-CD" => 0.6801, "GLH-CG" => -0.0174, "GLH-H" => 0.2719, "GLH-HA" => 0.0779, "GLH-HB2" => 0.0256, "GLH-HB3" => 0.0256, "GLH-HE2" => 0.4641, "GLH-HG2" => 0.043, "GLH-HG3" => 0.043, "GLH-N" => -0.4157, "GLH-O" => -0.5679, "GLH-OE1" => -0.5838, "GLH-OE2" => -0.6511,
        "GLN-C" => 0.5973, "GLN-CA" => -0.0031, "GLN-CB" => -0.0036, "GLN-CD" => 0.6951, "GLN-CG" => -0.0645, "GLN-H" => 0.2719, "GLN-HA" => 0.085, "GLN-HB2" => 0.0171, "GLN-HB3" => 0.0171, "GLN-HE21" => 0.4251, "GLN-HE22" => 0.4251, "GLN-HG2" => 0.0352, "GLN-HG3" => 0.0352, "GLN-N" => -0.4157, "GLN-NE2" => -0.9407, "GLN-O" => -0.5679, "GLN-OE1" => -0.6086,
        "GLU-C" => 0.5366, "GLU-CA" => 0.0397, "GLU-CB" => 0.056, "GLU-CD" => 0.8054, "GLU-CG" => 0.0136, "GLU-H" => 0.2936, "GLU-HA" => 0.1105, "GLU-HB2" => -0.0173, "GLU-HB3" => -0.0173, "GLU-HG2" => -0.0425, "GLU-HG3" => -0.0425, "GLU-N" => -0.5163, "GLU-O" => -0.5819, "GLU-OE1" => -0.8188, "GLU-OE2" => -0.8188,
        "GLY-C" => 0.5973, "GLY-CA" => -0.0252, "GLY-H" => 0.2719, "GLY-HA2" => 0.0698, "GLY-HA3" => 0.0698, "GLY-N" => -0.4157, "GLY-O" => -0.5679,
        "HID-C" => 0.5973, "HID-CA" => 0.0188, "HID-CB" => -0.0462, "HID-CD2" => 0.1292, "HID-CE1" => 0.2057, "HID-CG" => -0.0266, "HID-H" => 0.2719, "HID-HA" => 0.0881, "HID-HB2" => 0.0402, "HID-HB3" => 0.0402, "HID-HD1" => 0.3649, "HID-HD2" => 0.1147, "HID-HE1" => 0.1392, "HID-N" => -0.4157, "HID-ND1" => -0.3811, "HID-NE2" => -0.5727, "HID-O" => -0.5679,
        "HIS-C" => 0.5973, "HIS-CA" => 0.0188, "HIS-CB" => -0.0462, "HIS-CD2" => 0.1292, "HIS-CE1" => 0.2057, "HIS-CG" => -0.0266, "HIS-H" => 0.2719, "HIS-HA" => 0.0881, "HIS-HB2" => 0.0402, "HIS-HB3" => 0.0402, "HIS-HD1" => 0.3649, "HIS-HD2" => 0.1147, "HIS-HE1" => 0.1392, "HIS-HE2" => 0.3339, "HIS-N" => -0.4157, "HIS-ND1" => -0.3811, "HIS-NE2" => -0.5727, "HIS-O" => -0.5679,
        "HIE-C" => 0.5973, "HIE-CA" => -0.0581, "HIE-CB" => -0.0074, "HIE-CD2" => -0.2207, "HIE-CE1" => 0.1635, "HIE-CG" => 0.1868, "HIE-H" => 0.2719, "HIE-HA" => 0.136, "HIE-HB2" => 0.0367, "HIE-HB3" => 0.0367, "HIE-HD2" => 0.1862, "HIE-HE1" => 0.1435, "HIE-HE2" => 0.3339, "HIE-N" => -0.4157, "HIE-ND1" => -0.5432, "HIE-NE2" => -0.2795, "HIE-O" => -0.5679,
        "HIP-C" => 0.7341, "HIP-CA" => -0.1354, "HIP-CB" => -0.0414, "HIP-CD2" => -0.1141, "HIP-CE1" => -0.017, "HIP-CG" => -0.0012, "HIP-H" => 0.2747, "HIP-HA" => 0.1212, "HIP-HB2" => 0.081, "HIP-HB3" => 0.081, "HIP-HD1" => 0.3866, "HIP-HD2" => 0.2317, "HIP-HE1" => 0.2681, "HIP-HE2" => 0.3911, "HIP-N" => -0.3479, "HIP-ND1" => -0.1513, "HIP-NE2" => -0.1718, "HIP-O" => -0.5894,
        "ILE-C" => 0.5973, "ILE-CA" => -0.0597, "ILE-CB" => 0.1303, "ILE-CD1" => -0.066, "ILE-CG1" => -0.043, "ILE-CG2" => -0.3204, "ILE-H" => 0.2719, "ILE-HA" => 0.0869, "ILE-HB" => 0.0187, "ILE-HD11" => 0.0186, "ILE-HD12" => 0.0186, "ILE-HD13" => 0.0186, "ILE-HG12" => 0.0236, "ILE-HG13" => 0.0236, "ILE-HG21" => 0.0882, "ILE-HG22" => 0.0882, "ILE-HG23" => 0.0882, "ILE-N" => -0.4157, "ILE-O" => -0.5679,
        "LEU-C" => 0.5973, "LEU-CA" => -0.0518, "LEU-CB" => -0.1102, "LEU-CD1" => -0.4121, "LEU-CD2" => -0.4121, "LEU-CG" => 0.3531, "LEU-H" => 0.2719, "LEU-HA" => 0.0922, "LEU-HB2" => 0.0457, "LEU-HB3" => 0.0457, "LEU-HD11" => 0.1, "LEU-HD12" => 0.1, "LEU-HD13" => 0.1, "LEU-HD21" => 0.1, "LEU-HD22" => 0.1, "LEU-HD23" => 0.1, "LEU-HG" => -0.0361, "LEU-N" => -0.4157, "LEU-O" => -0.5679,
        "LYN-C" => 0.5973, "LYN-CA" => -0.07206, "LYN-CB" => -0.04845, "LYN-CD" => -0.03768, "LYN-CE" => 0.32604, "LYN-CG" => 0.06612, "LYN-H" => 0.2719, "LYN-HA" => 0.0994, "LYN-HB2" => 0.034, "LYN-HB3" => 0.034, "LYN-HD2" => 0.01155, "LYN-HD3" => 0.01155, "LYN-HE2" => -0.03358, "LYN-HE3" => -0.03358, "LYN-HG2" => 0.01041, "LYN-HG3" => 0.01041, "LYN-HZ2" => 0.38604, "LYN-HZ3" => 0.38604, "LYN-N" => -0.4157, "LYN-NZ" => -1.03581, "LYN-O" => -0.5679,
        "LYS-C" => 0.7341, "LYS-CA" => -0.24, "LYS-CB" => -0.0094, "LYS-CD" => -0.0479, "LYS-CE" => -0.0143, "LYS-CG" => 0.0187, "LYS-H" => 0.2747, "LYS-HA" => 0.1426, "LYS-HB2" => 0.0362, "LYS-HB3" => 0.0362, "LYS-HD2" => 0.0621, "LYS-HD3" => 0.0621, "LYS-HE2" => 0.1135, "LYS-HE3" => 0.1135, "LYS-HG2" => 0.0103, "LYS-HG3" => 0.0103, "LYS-HZ1" => 0.34, "LYS-HZ2" => 0.34, "LYS-HZ3" => 0.34, "LYS-N" => -0.3479, "LYS-NZ" => -0.3854, "LYS-O" => -0.5894,
        "MET-C" => 0.5973, "MET-CA" => -0.0237, "MET-CB" => 0.0342, "MET-CE" => -0.0536, "MET-CG" => 0.0018, "MET-H" => 0.2719, "MET-HA" => 0.088, "MET-HB2" => 0.0241, "MET-HB3" => 0.0241, "MET-HE1" => 0.0684, "MET-HE2" => 0.0684, "MET-HE3" => 0.0684, "MET-HG2" => 0.044, "MET-HG3" => 0.044, "MET-N" => -0.4157, "MET-O" => -0.5679, "MET-SD" => -0.2737,
        "PHE-C" => 0.5973, "PHE-CA" => -0.0024, "PHE-CB" => -0.0343, "PHE-CD1" => -0.1256, "PHE-CD2" => -0.1256, "PHE-CE1" => -0.1704, "PHE-CE2" => -0.1704, "PHE-CG" => 0.0118, "PHE-CZ" => -0.1072, "PHE-H" => 0.2719, "PHE-HA" => 0.0978, "PHE-HB2" => 0.0295, "PHE-HB3" => 0.0295, "PHE-HD1" => 0.133, "PHE-HD2" => 0.133, "PHE-HE1" => 0.143, "PHE-HE2" => 0.143, "PHE-HZ" => 0.1297, "PHE-N" => -0.4157, "PHE-O" => -0.5679,
        "PRO-C" => 0.5896, "PRO-CA" => -0.0266, "PRO-CB" => -0.007, "PRO-CD" => 0.0192, "PRO-CG" => 0.0189, "PRO-HA" => 0.0641, "PRO-HB2" => 0.0253, "PRO-HB3" => 0.0253, "PRO-HD2" => 0.0391, "PRO-HD3" => 0.0391, "PRO-HG2" => 0.0213, "PRO-HG3" => 0.0213, "PRO-N" => -0.2548, "PRO-O" => -0.5748,
        "RA-C1'" => 0.0394, "RA-C2" => 0.5875, "RA-C2'" => 0.067, "RA-C3'" => 0.2022, "RA-C4" => 0.3053, "RA-C4'" => 0.1065, "RA-C5" => 0.0515, "RA-C5'" => 0.0558, "RA-C6" => 0.7009, "RA-C8" => 0.2006, "RA-H1'" => 0.2007, "RA-H2" => 0.0473, "RA-H2'1" => 0.0972, "RA-H3'" => 0.0615, "RA-H4'" => 0.1174, "RA-H5'1" => 0.0679, "RA-H5'2" => 0.0679, "RA-H61" => 0.4115, "RA-H62" => 0.4115, "RA-H8" => 0.1553, "RA-HO'2" => 0.4186, "RA-N1" => -0.7615, "RA-N3" => -0.6997, "RA-N6" => -0.9019, "RA-N7" => -0.6073, "RA-N9" => -0.0251, "RA-O1P" => -0.776, "RA-O2'" => -0.6139, "RA-O2P" => -0.776, "RA-O3'" => -0.5246, "RA-O4'" => -0.3548, "RA-O5'" => -0.4989, "RA-P" => 1.1662,
        "RA3-C1'" => 0.0394, "RA3-C2" => 0.5875, "RA3-C2'" => 0.067, "RA3-C3'" => 0.2022, "RA3-C4" => 0.3053, "RA3-C4'" => 0.1065, "RA3-C5" => 0.0515, "RA3-C5'" => 0.0558, "RA3-C6" => 0.7009, "RA3-C8" => 0.2006, "RA3-H1'" => 0.2007, "RA3-H2" => 0.0473, "RA3-H2'1" => 0.0972, "RA3-H3'" => 0.0615, "RA3-H3T" => 0.4376, "RA3-H4'" => 0.1174, "RA3-H5'1" => 0.0679, "RA3-H5'2" => 0.0679, "RA3-H61" => 0.4115, "RA3-H62" => 0.4115, "RA3-H8" => 0.1553, "RA3-HO'2" => 0.4186, "RA3-N1" => -0.7615, "RA3-N3" => -0.6997, "RA3-N6" => -0.9019, "RA3-N7" => -0.6073, "RA3-N9" => -0.0251, "RA3-O1P" => -0.776, "RA3-O2'" => -0.6139, "RA3-O2P" => -0.776, "RA3-O3'" => -0.6541, "RA3-O4'" => -0.3548, "RA3-O5'" => -0.4989, "RA3-P" => 1.1662,
        "RA5-C1'" => 0.0394, "RA5-C2" => 0.5875, "RA5-C2'" => 0.067, "RA5-C3'" => 0.2022, "RA5-C4" => 0.3053, "RA5-C4'" => 0.1065, "RA5-C5" => 0.0515, "RA5-C5'" => 0.0558, "RA5-C6" => 0.7009, "RA5-C8" => 0.2006, "RA5-H1'" => 0.2007, "RA5-H2" => 0.0473, "RA5-H2'1" => 0.0972, "RA5-H3'" => 0.0615, "RA5-H4'" => 0.1174, "RA5-H5'1" => 0.0679, "RA5-H5'2" => 0.0679, "RA5-H5T" => 0.4295, "RA5-H61" => 0.4115, "RA5-H62" => 0.4115, "RA5-H8" => 0.1553, "RA5-HO'2" => 0.4186, "RA5-N1" => -0.7615, "RA5-N3" => -0.6997, "RA5-N6" => -0.9019, "RA5-N7" => -0.6073, "RA5-N9" => -0.0251, "RA5-O2'" => -0.6139, "RA5-O3'" => -0.5246, "RA5-O4'" => -0.3548, "RA5-O5'" => -0.6223,
        "RAN-C1'" => 0.0394, "RAN-C2" => 0.5875, "RAN-C2'" => 0.067, "RAN-C3'" => 0.2022, "RAN-C4" => 0.3053, "RAN-C4'" => 0.1065, "RAN-C5" => 0.0515, "RAN-C5'" => 0.0558, "RAN-C6" => 0.7009, "RAN-C8" => 0.2006, "RAN-H1'" => 0.2007, "RAN-H2" => 0.0473, "RAN-H2'1" => 0.0972, "RAN-H3'" => 0.0615, "RAN-H3T" => 0.4376, "RAN-H4'" => 0.1174, "RAN-H5'1" => 0.0679, "RAN-H5'2" => 0.0679, "RAN-H5T" => 0.4295, "RAN-H61" => 0.4115, "RAN-H62" => 0.4115, "RAN-H8" => 0.1553, "RAN-HO'2" => 0.4186, "RAN-N1" => -0.7615, "RAN-N3" => -0.6997, "RAN-N6" => -0.9019, "RAN-N7" => -0.6073, "RAN-N9" => -0.0251, "RAN-O2'" => -0.6139, "RAN-O3'" => -0.6541, "RAN-O4'" => -0.3548, "RAN-O5'" => -0.6223,
        "RC-C1'" => 0.0066, "RC-C2" => 0.7538, "RC-C2'" => 0.067, "RC-C3'" => 0.2022, "RC-C4" => 0.8185, "RC-C4'" => 0.1065, "RC-C5" => -0.5215, "RC-C5'" => 0.0558, "RC-C6" => 0.0053, "RC-H1'" => 0.2029, "RC-H2'1" => 0.0972, "RC-H3'" => 0.0615, "RC-H4'" => 0.1174, "RC-H41" => 0.4234, "RC-H42" => 0.4234, "RC-H5" => 0.1928, "RC-H5'1" => 0.0679, "RC-H5'2" => 0.0679, "RC-H6" => 0.1958, "RC-HO'2" => 0.4186, "RC-N1" => -0.0484, "RC-N3" => -0.7584, "RC-N4" => -0.953, "RC-O1P" => -0.776, "RC-O2" => -0.6252, "RC-O2'" => -0.6139, "RC-O2P" => -0.776, "RC-O3'" => -0.5246, "RC-O4'" => -0.3548, "RC-O5'" => -0.4989, "RC-P" => 1.1662,
        "RC3-C1'" => 0.0066, "RC3-C2" => 0.7538, "RC3-C2'" => 0.067, "RC3-C3'" => 0.2022, "RC3-C4" => 0.8185, "RC3-C4'" => 0.1065, "RC3-C5" => -0.5215, "RC3-C5'" => 0.0558, "RC3-C6" => 0.0053, "RC3-H1'" => 0.2029, "RC3-H2'1" => 0.0972, "RC3-H3'" => 0.0615, "RC3-H3T" => 0.4376, "RC3-H4'" => 0.1174, "RC3-H41" => 0.4234, "RC3-H42" => 0.4234, "RC3-H5" => 0.1928, "RC3-H5'1" => 0.0679, "RC3-H5'2" => 0.0679, "RC3-H6" => 0.1958, "RC3-HO'2" => 0.4186, "RC3-N1" => -0.0484, "RC3-N3" => -0.7584, "RC3-N4" => -0.953, "RC3-O1P" => -0.776, "RC3-O2" => -0.6252, "RC3-O2'" => -0.6139, "RC3-O2P" => -0.776, "RC3-O3'" => -0.6541, "RC3-O4'" => -0.3548, "RC3-O5'" => -0.4989, "RC3-P" => 1.1662,
        "RC5-C1'" => 0.0066, "RC5-C2" => 0.7538, "RC5-C2'" => 0.067, "RC5-C3'" => 0.2022, "RC5-C4" => 0.8185, "RC5-C4'" => 0.1065, "RC5-C5" => -0.5215, "RC5-C5'" => 0.0558, "RC5-C6" => 0.0053, "RC5-H1'" => 0.2029, "RC5-H2'1" => 0.0972, "RC5-H3'" => 0.0615, "RC5-H4'" => 0.1174, "RC5-H41" => 0.4234, "RC5-H42" => 0.4234, "RC5-H5" => 0.1928, "RC5-H5'1" => 0.0679, "RC5-H5'2" => 0.0679, "RC5-H5T" => 0.4295, "RC5-H6" => 0.1958, "RC5-HO'2" => 0.4186, "RC5-N1" => -0.0484, "RC5-N3" => -0.7584, "RC5-N4" => -0.953, "RC5-O2" => -0.6252, "RC5-O2'" => -0.6139, "RC5-O3'" => -0.5246, "RC5-O4'" => -0.3548, "RC5-O5'" => -0.6223,
        "RCN-C1'" => 0.0066, "RCN-C2" => 0.7538, "RCN-C2'" => 0.067, "RCN-C3'" => 0.2022, "RCN-C4" => 0.8185, "RCN-C4'" => 0.1065, "RCN-C5" => -0.5215, "RCN-C5'" => 0.0558, "RCN-C6" => 0.0053, "RCN-H1'" => 0.2029, "RCN-H2'1" => 0.0972, "RCN-H3'" => 0.0615, "RCN-H3T" => 0.4376, "RCN-H4'" => 0.1174, "RCN-H41" => 0.4234, "RCN-H42" => 0.4234, "RCN-H5" => 0.1928, "RCN-H5'1" => 0.0679, "RCN-H5'2" => 0.0679, "RCN-H5T" => 0.4295, "RCN-H6" => 0.1958, "RCN-HO'2" => 0.4186, "RCN-N1" => -0.0484, "RCN-N3" => -0.7584, "RCN-N4" => -0.953, "RCN-O2" => -0.6252, "RCN-O2'" => -0.6139, "RCN-O3'" => -0.6541, "RCN-O4'" => -0.3548, "RCN-O5'" => -0.6223,
        "RG-C1'" => 0.0191, "RG-C2" => 0.7657, "RG-C2'" => 0.067, "RG-C3'" => 0.2022, "RG-C4" => 0.1222, "RG-C4'" => 0.1065, "RG-C5" => 0.1744, "RG-C5'" => 0.0558, "RG-C6" => 0.477, "RG-C8" => 0.1374, "RG-H1" => 0.3424, "RG-H1'" => 0.2006, "RG-H2'1" => 0.0972, "RG-H21" => 0.4364, "RG-H22" => 0.4364, "RG-H3'" => 0.0615, "RG-H4'" => 0.1174, "RG-H5'1" => 0.0679, "RG-H5'2" => 0.0679, "RG-H8" => 0.164, "RG-HO'2" => 0.4186, "RG-N1" => -0.4787, "RG-N2" => -0.9672, "RG-N3" => -0.6323, "RG-N7" => -0.5709, "RG-N9" => 0.0492, "RG-O1P" => -0.776, "RG-O2'" => -0.6139, "RG-O2P" => -0.776, "RG-O3'" => -0.5246, "RG-O4'" => -0.3548, "RG-O5'" => -0.4989, "RG-O6" => -0.5597, "RG-P" => 1.1662,
        "RG3-C1'" => 0.0191, "RG3-C2" => 0.7657, "RG3-C2'" => 0.067, "RG3-C3'" => 0.2022, "RG3-C4" => 0.1222, "RG3-C4'" => 0.1065, "RG3-C5" => 0.1744, "RG3-C5'" => 0.0558, "RG3-C6" => 0.477, "RG3-C8" => 0.1374, "RG3-H1" => 0.3424, "RG3-H1'" => 0.2006, "RG3-H2'1" => 0.0972, "RG3-H21" => 0.4364, "RG3-H22" => 0.4364, "RG3-H3'" => 0.0615, "RG3-H3T" => 0.4376, "RG3-H4'" => 0.1174, "RG3-H5'1" => 0.0679, "RG3-H5'2" => 0.0679, "RG3-H8" => 0.164, "RG3-HO'2" => 0.4186, "RG3-N1" => -0.4787, "RG3-N2" => -0.9672, "RG3-N3" => -0.6323, "RG3-N7" => -0.5709, "RG3-N9" => 0.0492, "RG3-O1P" => -0.776, "RG3-O2'" => -0.6139, "RG3-O2P" => -0.776, "RG3-O3'" => -0.6541, "RG3-O4'" => -0.3548, "RG3-O5'" => -0.4989, "RG3-O6" => -0.5597, "RG3-P" => 1.1662,
        "RG5-C1'" => 0.0191, "RG5-C2" => 0.7657, "RG5-C2'" => 0.067, "RG5-C3'" => 0.2022, "RG5-C4" => 0.1222, "RG5-C4'" => 0.1065, "RG5-C5" => 0.1744, "RG5-C5'" => 0.0558, "RG5-C6" => 0.477, "RG5-C8" => 0.1374, "RG5-H1" => 0.3424, "RG5-H1'" => 0.2006, "RG5-H2'1" => 0.0972, "RG5-H21" => 0.4364, "RG5-H22" => 0.4364, "RG5-H3'" => 0.0615, "RG5-H4'" => 0.1174, "RG5-H5'1" => 0.0679, "RG5-H5'2" => 0.0679, "RG5-H5T" => 0.4295, "RG5-H8" => 0.164, "RG5-HO'2" => 0.4186, "RG5-N1" => -0.4787, "RG5-N2" => -0.9672, "RG5-N3" => -0.6323, "RG5-N7" => -0.5709, "RG5-N9" => 0.0492, "RG5-O2'" => -0.6139, "RG5-O3'" => -0.5246, "RG5-O4'" => -0.3548, "RG5-O5'" => -0.6223, "RG5-O6" => -0.5597,
        "RGN-C1'" => 0.0191, "RGN-C2" => 0.7657, "RGN-C2'" => 0.067, "RGN-C3'" => 0.2022, "RGN-C4" => 0.1222, "RGN-C4'" => 0.1065, "RGN-C5" => 0.1744, "RGN-C5'" => 0.0558, "RGN-C6" => 0.477, "RGN-C8" => 0.1374, "RGN-H1" => 0.3424, "RGN-H1'" => 0.2006, "RGN-H2'1" => 0.0972, "RGN-H21" => 0.4364, "RGN-H22" => 0.4364, "RGN-H3'" => 0.0615, "RGN-H3T" => 0.4376, "RGN-H4'" => 0.1174, "RGN-H5'1" => 0.0679, "RGN-H5'2" => 0.0679, "RGN-H5T" => 0.4295, "RGN-H8" => 0.164, "RGN-HO'2" => 0.4186, "RGN-N1" => -0.4787, "RGN-N2" => -0.9672, "RGN-N3" => -0.6323, "RGN-N7" => -0.5709, "RGN-N9" => 0.0492, "RGN-O2'" => -0.6139, "RGN-O3'" => -0.6541, "RGN-O4'" => -0.3548, "RGN-O5'" => -0.6223, "RGN-O6" => -0.5597,
        "RU-C1'" => 0.0674, "RU-C2" => 0.4687, "RU-C2'" => 0.067, "RU-C3'" => 0.2022, "RU-C4" => 0.5952, "RU-C4'" => 0.1065, "RU-C5" => -0.3635, "RU-C5'" => 0.0558, "RU-C6" => -0.1126, "RU-H1'" => 0.1824, "RU-H2'1" => 0.0972, "RU-H3" => 0.3154, "RU-H3'" => 0.0615, "RU-H4'" => 0.1174, "RU-H5" => 0.1811, "RU-H5'1" => 0.0679, "RU-H5'2" => 0.0679, "RU-H6" => 0.2188, "RU-HO'2" => 0.4186, "RU-N1" => 0.0418, "RU-N3" => -0.3549, "RU-O1P" => -0.776, "RU-O2" => -0.5477, "RU-O2'" => -0.6139, "RU-O2P" => -0.776, "RU-O3'" => -0.5246, "RU-O4" => -0.5761, "RU-O4'" => -0.3548, "RU-O5'" => -0.4989, "RU-P" => 1.1662,
        "RU3-C1'" => 0.0674, "RU3-C2" => 0.4687, "RU3-C2'" => 0.067, "RU3-C3'" => 0.2022, "RU3-C4" => 0.5952, "RU3-C4'" => 0.1065, "RU3-C5" => -0.3635, "RU3-C5'" => 0.0558, "RU3-C6" => -0.1126, "RU3-H1'" => 0.1824, "RU3-H2'1" => 0.0972, "RU3-H3" => 0.3154, "RU3-H3'" => 0.0615, "RU3-H3T" => 0.4376, "RU3-H4'" => 0.1174, "RU3-H5" => 0.1811, "RU3-H5'1" => 0.0679, "RU3-H5'2" => 0.0679, "RU3-H6" => 0.2188, "RU3-HO'2" => 0.4186, "RU3-N1" => 0.0418, "RU3-N3" => -0.3549, "RU3-O1P" => -0.776, "RU3-O2" => -0.5477, "RU3-O2'" => -0.6139, "RU3-O2P" => -0.776, "RU3-O3'" => -0.6541, "RU3-O4" => -0.5761, "RU3-O4'" => -0.3548, "RU3-O5'" => -0.4989, "RU3-P" => 1.1662,
        "RU5-C1'" => 0.0674, "RU5-C2" => 0.4687, "RU5-C2'" => 0.067, "RU5-C3'" => 0.2022, "RU5-C4" => 0.5952, "RU5-C4'" => 0.1065, "RU5-C5" => -0.3635, "RU5-C5'" => 0.0558, "RU5-C6" => -0.1126, "RU5-H1'" => 0.1824, "RU5-H2'1" => 0.0972, "RU5-H3" => 0.3154, "RU5-H3'" => 0.0615, "RU5-H4'" => 0.1174, "RU5-H5" => 0.1811, "RU5-H5'1" => 0.0679, "RU5-H5'2" => 0.0679, "RU5-H5T" => 0.4295, "RU5-H6" => 0.2188, "RU5-HO'2" => 0.4186, "RU5-N1" => 0.0418, "RU5-N3" => -0.3549, "RU5-O2" => -0.5477, "RU5-O2'" => -0.6139, "RU5-O3'" => -0.5246, "RU5-O4" => -0.5761, "RU5-O4'" => -0.3548, "RU5-O5'" => -0.6223,
        "RUN-C1'" => 0.0674, "RUN-C2" => 0.4687, "RUN-C2'" => 0.067, "RUN-C3'" => 0.2022, "RUN-C4" => 0.5952, "RUN-C4'" => 0.1065, "RUN-C5" => -0.3635, "RUN-C5'" => 0.0558, "RUN-C6" => -0.1126, "RUN-H1'" => 0.1824, "RUN-H2'1" => 0.0972, "RUN-H3" => 0.3154, "RUN-H3'" => 0.0615, "RUN-H3T" => 0.4376, "RUN-H4'" => 0.1174, "RUN-H5" => 0.1811, "RUN-H5'1" => 0.0679, "RUN-H5'2" => 0.0679, "RUN-H5T" => 0.4295, "RUN-H6" => 0.2188, "RUN-HO'2" => 0.4186, "RUN-N1" => 0.0418, "RUN-N3" => -0.3549, "RUN-O2" => -0.5477, "RUN-O2'" => -0.6139, "RUN-O3'" => -0.6541, "RUN-O4" => -0.5761, "RUN-O4'" => -0.3548, "RUN-O5'" => -0.6223,
        "SER-C" => 0.5973, "SER-CA" => -0.0249, "SER-CB" => 0.2117, "SER-H" => 0.2719, "SER-HA" => 0.0843, "SER-HB2" => 0.0352, "SER-HB3" => 0.0352, "SER-HG" => 0.4275, "SER-N" => -0.4157, "SER-O" => -0.5679, "SER-OG" => -0.6546,
        "THR-C" => 0.5973, "THR-CA" => -0.0389, "THR-CB" => 0.3654, "THR-CG2" => -0.2438, "THR-H" => 0.2719, "THR-HA" => 0.1007, "THR-HB" => 0.0043, "THR-HG1" => 0.4102, "THR-HG21" => 0.0642, "THR-HG22" => 0.0642, "THR-HG23" => 0.0642, "THR-N" => -0.4157, "THR-O" => -0.5679, "THR-OG1" => -0.6761,
        "TRP-C" => 0.5973, "TRP-CA" => -0.0275, "TRP-CB" => -0.005, "TRP-CD1" => -0.1638, "TRP-CD2" => 0.1243, "TRP-CE2" => 0.138, "TRP-CE3" => -0.2387, "TRP-CG" => -0.1415, "TRP-CH2" => -0.1134, "TRP-CZ2" => -0.2601, "TRP-CZ3" => -0.1972, "TRP-H" => 0.2719, "TRP-HA" => 0.1123, "TRP-HB2" => 0.0339, "TRP-HB3" => 0.0339, "TRP-HD1" => 0.2062, "TRP-HE1" => 0.3412, "TRP-HE3" => 0.17, "TRP-HH2" => 0.1417, "TRP-HZ2" => 0.1572, "TRP-HZ3" => 0.1447, "TRP-N" => -0.4157, "TRP-NE1" => -0.3418, "TRP-O" => -0.5679,
        "TYR-C" => 0.5973, "TYR-CA" => -0.0014, "TYR-CB" => -0.0152, "TYR-CD1" => -0.1906, "TYR-CD2" => -0.1906, "TYR-CE1" => -0.2341, "TYR-CE2" => -0.2341, "TYR-CG" => -0.0011, "TYR-CZ" => 0.3226, "TYR-H" => 0.2719, "TYR-HA" => 0.0876, "TYR-HB2" => 0.0295, "TYR-HB3" => 0.0295, "TYR-HD1" => 0.1699, "TYR-HD2" => 0.1699, "TYR-HE1" => 0.1656, "TYR-HE2" => 0.1656, "TYR-HH" => 0.3992, "TYR-N" => -0.4157, "TYR-O" => -0.5679, "TYR-OH" => -0.5579,
        "VAL-C" => 0.5973, "VAL-CA" => -0.0875, "VAL-CB" => 0.2985, "VAL-CG1" => -0.3192, "VAL-CG2" => -0.3192, "VAL-H" => 0.2719, "VAL-HA" => 0.0969, "VAL-HB" => -0.0297, "VAL-HG11" => 0.0791, "VAL-HG12" => 0.0791, "VAL-HG13" => 0.0791, "VAL-HG21" => 0.0791, "VAL-HG22" => 0.0791, "VAL-HG23" => 0.0791, "VAL-N" => -0.4157, "VAL-O" => -0.5679];

    static ref NT_ELE_CHARGES: HashMap<&'static str, f64> = hashmap![
        "ACE-C" => 0.5972, "ACE-CH3" => -0.3662, "ACE-HH31" => 0.1123, "ACE-HH32" => 0.1123, "ACE-HH33" => 0.1123, "ACE-O" => -0.5679,
        "ALA-C" => 0.6163, "ALA-CA" => 0.0962, "ALA-CB" => -0.0597, "ALA-H1" => 0.1997, "ALA-H2" => 0.1997, "ALA-H3" => 0.1997, "ALA-HA" => 0.0889, "ALA-HB1" => 0.03, "ALA-HB2" => 0.03, "ALA-HB3" => 0.03, "ALA-N" => 0.1414, "ALA-O" => -0.5722,
        "ARG-C" => 0.7214, "ARG-CA" => -0.0223, "ARG-CB" => 0.0118, "ARG-CD" => 0.0935, "ARG-CG" => 0.0236, "ARG-CZ" => 0.8281, "ARG-H1" => 0.2083, "ARG-H2" => 0.2083, "ARG-H3" => 0.2083, "ARG-HA" => 0.1242, "ARG-HB2" => 0.0226, "ARG-HB3" => 0.0226, "ARG-HD2" => 0.0527, "ARG-HD3" => 0.0527, "ARG-HE" => 0.3592, "ARG-HG2" => 0.0309, "ARG-HG3" => 0.0309, "ARG-HH11" => 0.4494, "ARG-HH12" => 0.4494, "ARG-HH21" => 0.4494, "ARG-HH22" => 0.4494, "ARG-N" => 0.1305, "ARG-NE" => -0.565, "ARG-NH1" => -0.8693, "ARG-NH2" => -0.8693, "ARG-O" => -0.6013,
        "ASN-C" => 0.6163, "ASN-CA" => 0.0368, "ASN-CB" => -0.0283, "ASN-CG" => 0.5833, "ASN-H1" => 0.1921, "ASN-H2" => 0.1921, "ASN-H3" => 0.1921, "ASN-HA" => 0.1231, "ASN-HB2" => 0.0515, "ASN-HB3" => 0.0515, "ASN-HD21" => 0.4097, "ASN-HD22" => 0.4097, "ASN-N" => 0.1801, "ASN-ND2" => -0.8634, "ASN-O" => -0.5722, "ASN-OD1" => -0.5744,
        "ASP-C" => 0.5621, "ASP-CA" => 0.0292, "ASP-CB" => -0.0235, "ASP-CG" => 0.8194, "ASP-H1" => 0.22, "ASP-H2" => 0.22, "ASP-H3" => 0.22, "ASP-HA" => 0.1141, "ASP-HB2" => -0.0169, "ASP-HB3" => -0.0169, "ASP-N" => 0.0782, "ASP-O" => -0.5889, "ASP-OD1" => -0.8084, "ASP-OD2" => -0.8084,
        "CYS-C" => 0.6123, "CYS-CA" => 0.0927, "CYS-CB" => -0.1195, "CYS-H1" => 0.2023, "CYS-H2" => 0.2023, "CYS-H3" => 0.2023, "CYS-HA" => 0.1411, "CYS-HB2" => 0.1188, "CYS-HB3" => 0.1188, "CYS-HSG" => 0.1975, "CYS-N" => 0.1325, "CYS-O" => -0.5713, "CYS-SG" => -0.3298,
        "CYX-C" => 0.6123, "CYX-CA" => 0.1055, "CYX-CB" => -0.0277, "CYX-H1" => 0.1815, "CYX-H2" => 0.1815, "CYX-H3" => 0.1815, "CYX-HA" => 0.0922, "CYX-HB2" => 0.068, "CYX-HB3" => 0.068, "CYX-N" => 0.2069, "CYX-O" => -0.5713, "CYX-SG" => -0.0984,
        "GLN-C" => 0.6123, "GLN-CA" => 0.0536, "GLN-CB" => 0.0651, "GLN-CD" => 0.7354, "GLN-CG" => -0.0903, "GLN-H1" => 0.1996, "GLN-H2" => 0.1996, "GLN-H3" => 0.1996, "GLN-HA" => 0.1015, "GLN-HB2" => 0.005, "GLN-HB3" => 0.005, "GLN-HE21" => 0.4429, "GLN-HE22" => 0.4429, "GLN-HG2" => 0.0331, "GLN-HG3" => 0.0331, "GLN-N" => 0.1493, "GLN-NE2" => -1.0031, "GLN-O" => -0.5713, "GLN-OE1" => -0.6133,
        "GLU-C" => 0.5621, "GLU-CA" => 0.0588, "GLU-CB" => 0.0909, "GLU-CD" => 0.8087, "GLU-CG" => -0.0236, "GLU-H1" => 0.2391, "GLU-H2" => 0.2391, "GLU-H3" => 0.2391, "GLU-HA" => 0.1202, "GLU-HB2" => -0.0232, "GLU-HB3" => -0.0232, "GLU-HG2" => -0.0315, "GLU-HG3" => -0.0315, "GLU-N" => 0.0017, "GLU-O" => -0.5889, "GLU-OE1" => -0.8189, "GLU-OE2" => -0.8189,
        "GLY-C" => 0.6163, "GLY-CA" => -0.01, "GLY-H1" => 0.1642, "GLY-H2" => 0.1642, "GLY-H3" => 0.1642, "GLY-HA2" => 0.0895, "GLY-HA3" => 0.0895, "GLY-N" => 0.2943, "GLY-O" => -0.5722,
        "HID-C" => 0.6123, "HID-CA" => 0.0964, "HID-CB" => 0.0259, "HID-CD2" => 0.1046, "HID-CE1" => 0.2127, "HID-CG" => -0.0399, "HID-H1" => 0.1963, "HID-H2" => 0.1963, "HID-H3" => 0.1963, "HID-HA" => 0.0958, "HID-HB2" => 0.0209, "HID-HB3" => 0.0209, "HID-HD1" => 0.3632, "HID-HD2" => 0.1299, "HID-HE1" => 0.1385, "HID-N" => 0.1542, "HID-ND1" => -0.3819, "HID-NE2" => -0.5711, "HID-O" => -0.5713,
        "HIS-C" => 0.6123, "HIS-CA" => 0.0964, "HIS-CB" => 0.0259, "HIS-CD2" => 0.1046, "HIS-CE1" => 0.2127, "HIS-CG" => -0.0399, "HIS-H1" => 0.1963, "HIS-H2" => 0.1963, "HIS-H3" => 0.1963, "HIS-HA" => 0.0958, "HIS-HB2" => 0.0209, "HIS-HB3" => 0.0209, "HIS-HD1" => 0.3632, "HIS-HD2" => 0.1299, "HIS-HE1" => 0.1385, "HIS-HE2" => 0.3324, "HIS-N" => 0.1542, "HIS-ND1" => -0.3819, "HIS-NE2" => -0.5711, "HIS-O" => -0.5713,
        "HIE-C" => 0.6123, "HIE-CA" => 0.0236, "HIE-CB" => 0.0489, "HIE-CD2" => -0.2349, "HIE-CE1" => 0.1804, "HIE-CG" => 0.174, "HIE-H1" => 0.2016, "HIE-H2" => 0.2016, "HIE-H3" => 0.2016, "HIE-HA" => 0.138, "HIE-HB2" => 0.0223, "HIE-HB3" => 0.0223, "HIE-HD2" => 0.1963, "HIE-HE1" => 0.1397, "HIE-HE2" => 0.3324, "HIE-N" => 0.1472, "HIE-ND1" => -0.5579, "HIE-NE2" => -0.2781, "HIE-O" => -0.5713,
        "HIP-C" => 0.7214, "HIP-CA" => 0.0581, "HIP-CB" => 0.0484, "HIP-CD2" => -0.1433, "HIP-CE1" => -0.0011, "HIP-CG" => -0.0236, "HIP-H1" => 0.1704, "HIP-H2" => 0.1704, "HIP-H3" => 0.1704, "HIP-HA" => 0.1047, "HIP-HB2" => 0.0531, "HIP-HB3" => 0.0531, "HIP-HD1" => 0.3821, "HIP-HD2" => 0.2495, "HIP-HE1" => 0.2645, "HIP-HE2" => 0.3921, "HIP-N" => 0.256, "HIP-ND1" => -0.151, "HIP-NE2" => -0.1739, "HIP-O" => -0.6013,
        "ILE-C" => 0.6123, "ILE-CA" => 0.0257, "ILE-CB" => 0.1885, "ILE-CD1" => -0.0908, "ILE-CG1" => -0.0387, "ILE-CG2" => -0.372, "ILE-H1" => 0.2329, "ILE-H2" => 0.2329, "ILE-H3" => 0.2329, "ILE-HA" => 0.1031, "ILE-HB" => 0.0213, "ILE-HD11" => 0.0226, "ILE-HD12" => 0.0226, "ILE-HD13" => 0.0226, "ILE-HG12" => 0.0201, "ILE-HG13" => 0.0201, "ILE-HG21" => 0.0947, "ILE-HG22" => 0.0947, "ILE-HG23" => 0.0947, "ILE-N" => 0.0311, "ILE-O" => -0.5713,
        "LEU-C" => 0.6123, "LEU-CA" => 0.0104, "LEU-CB" => -0.0244, "LEU-CD1" => -0.4106, "LEU-CD2" => -0.4104, "LEU-CG" => 0.3421, "LEU-H1" => 0.2148, "LEU-H2" => 0.2148, "LEU-H3" => 0.2148, "LEU-HA" => 0.1053, "LEU-HB2" => 0.0256, "LEU-HB3" => 0.0256, "LEU-HD11" => 0.098, "LEU-HD12" => 0.098, "LEU-HD13" => 0.098, "LEU-HD21" => 0.098, "LEU-HD22" => 0.098, "LEU-HD23" => 0.098, "LEU-HG" => -0.038, "LEU-N" => 0.101, "LEU-O" => -0.5713,
        "LYS-C" => 0.7214, "LYS-CA" => -0.0015, "LYS-CB" => 0.0212, "LYS-CD" => -0.0608, "LYS-CE" => -0.0181, "LYS-CG" => -0.0048, "LYS-H1" => 0.2165, "LYS-H2" => 0.2165, "LYS-H3" => 0.2165, "LYS-HA" => 0.118, "LYS-HB2" => 0.0283, "LYS-HB3" => 0.0283, "LYS-HD2" => 0.0633, "LYS-HD3" => 0.0633, "LYS-HE2" => 0.1171, "LYS-HE3" => 0.1171, "LYS-HG2" => 0.0121, "LYS-HG3" => 0.0121, "LYS-HZ1" => 0.3382, "LYS-HZ2" => 0.3382, "LYS-HZ3" => 0.3382, "LYS-N" => 0.0966, "LYS-NZ" => -0.3764, "LYS-O" => -0.6013,
        "MET-C" => 0.6123, "MET-CA" => 0.0221, "MET-CB" => 0.0865, "MET-CE" => -0.0341, "MET-CG" => 0.0334, "MET-H1" => 0.1984, "MET-H2" => 0.1984, "MET-H3" => 0.1984, "MET-HA" => 0.1116, "MET-HB2" => 0.0125, "MET-HB3" => 0.0125, "MET-HE1" => 0.0597, "MET-HE2" => 0.0597, "MET-HE3" => 0.0597, "MET-HG2" => 0.0292, "MET-HG3" => 0.0292, "MET-N" => 0.1592, "MET-O" => -0.5713, "MET-SD" => -0.2774,
        "PHE-C" => 0.6123, "PHE-CA" => 0.0733, "PHE-CB" => 0.033, "PHE-CD1" => -0.1392, "PHE-CD2" => -0.1391, "PHE-CE1" => -0.1602, "PHE-CE2" => -0.1603, "PHE-CG" => 0.0031, "PHE-CZ" => -0.1208, "PHE-H1" => 0.1921, "PHE-H2" => 0.1921, "PHE-H3" => 0.1921, "PHE-HA" => 0.1041, "PHE-HB2" => 0.0104, "PHE-HB3" => 0.0104, "PHE-HD1" => 0.1374, "PHE-HD2" => 0.1374, "PHE-HE1" => 0.1433, "PHE-HE2" => 0.1433, "PHE-HZ" => 0.1329, "PHE-N" => 0.1737, "PHE-O" => -0.5713,
        "PRO-C" => 0.526, "PRO-CA" => 0.1, "PRO-CB" => -0.115, "PRO-CD" => -0.012, "PRO-CG" => -0.121, "PRO-H2" => 0.312, "PRO-H3" => 0.312, "PRO-HA" => 0.1, "PRO-HB2" => 0.1, "PRO-HB3" => 0.1, "PRO-HD2" => 0.1, "PRO-HD3" => 0.1, "PRO-HG2" => 0.1, "PRO-HG3" => 0.1, "PRO-N" => -0.202, "PRO-O" => -0.5,
        "SER-C" => 0.6163, "SER-CA" => 0.0567, "SER-CB" => 0.2596, "SER-H1" => 0.1898, "SER-H2" => 0.1898, "SER-H3" => 0.1898, "SER-HA" => 0.0782, "SER-HB2" => 0.0273, "SER-HB3" => 0.0273, "SER-HG" => 0.4239, "SER-N" => 0.1849, "SER-O" => -0.5722, "SER-OG" => -0.6714,
        "THR-C" => 0.6163, "THR-CA" => 0.0034, "THR-CB" => 0.4514, "THR-CG2" => -0.2554, "THR-H1" => 0.1934, "THR-H2" => 0.1934, "THR-H3" => 0.1934, "THR-HA" => 0.1087, "THR-HB" => -0.0323, "THR-HG1" => 0.407, "THR-HG21" => 0.0627, "THR-HG22" => 0.0627, "THR-HG23" => 0.0627, "THR-N" => 0.1812, "THR-O" => -0.5722, "THR-OG1" => -0.6764,
        "TRP-C" => 0.6123, "TRP-CA" => 0.0421, "TRP-CB" => 0.0543, "TRP-CD1" => -0.1788, "TRP-CD2" => 0.1132, "TRP-CE2" => 0.1575, "TRP-CE3" => -0.2265, "TRP-CG" => -0.1654, "TRP-CH2" => -0.108, "TRP-CZ2" => -0.271, "TRP-CZ3" => -0.2034, "TRP-H1" => 0.1888, "TRP-H2" => 0.1888, "TRP-H3" => 0.1888, "TRP-HA" => 0.1162, "TRP-HB2" => 0.0222, "TRP-HB3" => 0.0222, "TRP-HD1" => 0.2195, "TRP-HE1" => 0.3412, "TRP-HE3" => 0.1646, "TRP-HH2" => 0.1411, "TRP-HZ2" => 0.1589, "TRP-HZ3" => 0.1458, "TRP-N" => 0.1913, "TRP-NE1" => -0.3444, "TRP-O" => -0.5713,
        "TYR-C" => 0.6123, "TYR-CA" => 0.057, "TYR-CB" => 0.0659, "TYR-CD1" => -0.2002, "TYR-CD2" => -0.2002, "TYR-CE1" => -0.2239, "TYR-CE2" => -0.2239, "TYR-CG" => -0.0205, "TYR-CZ" => 0.3139, "TYR-H1" => 0.1873, "TYR-H2" => 0.1873, "TYR-H3" => 0.1873, "TYR-HA" => 0.0983, "TYR-HB2" => 0.0102, "TYR-HB3" => 0.0102, "TYR-HD1" => 0.172, "TYR-HD2" => 0.172, "TYR-HE1" => 0.165, "TYR-HE2" => 0.165, "TYR-HH" => 0.4001, "TYR-N" => 0.194, "TYR-O" => -0.5713, "TYR-OH" => -0.5578,
        "VAL-C" => 0.6163, "VAL-CA" => -0.0054, "VAL-CB" => 0.3196, "VAL-CG1" => -0.3129, "VAL-CG2" => -0.3129, "VAL-H1" => 0.2272, "VAL-H2" => 0.2272, "VAL-H3" => 0.2272, "VAL-HA" => 0.1093, "VAL-HB" => -0.0221, "VAL-HG11" => 0.0735, "VAL-HG12" => 0.0735, "VAL-HG13" => 0.0735, "VAL-HG21" => 0.0735, "VAL-HG22" => 0.0735, "VAL-HG23" => 0.0735, "VAL-N" => 0.0577, "VAL-O" => -0.5722];
}

pub struct DNADockingModel {
    pub atoms: Vec<usize>,
    pub coordinates: Vec<[f64; 3]>,
    pub membrane: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub num_anm: usize,
    pub nmodes: Vec<f64>,
    pub vdw_radii: Vec<f64>,
    pub vdw_charges: Vec<f64>,
    pub ele_charges: Vec<f64>,
}

impl<'a> DNADockingModel {

    fn new(structure: &'a Structure, active_restraints: &'a [String], passive_restraints: &'a [String],
        nmodes: &[f64], num_anm: usize) -> DNADockingModel {
        let mut model = DNADockingModel {
            atoms: Vec::new(),
            coordinates: Vec::new(),
            membrane: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
            nmodes: nmodes.to_owned(),
            num_anm,
            vdw_radii: Vec::new(),
            vdw_charges: Vec::new(),
            ele_charges: Vec::new(),
        };

        let mut atom_index: u64 = 0;
        for chain in structure.chains.iter() {
            for residue in chain.lst_res.iter() {
                let mut res_id = format!("{}.{}.{}", chain.name, residue.name.trim(), residue.res_num);
                if let Some(c) = residue.res_icode {
                    res_id.push(c);
                }

                for atom in residue.lst_atom.iter() {
                    // Membrane beads MMB.BJ
                    let rec_atom_type = format!("{}{}", residue.name.trim(), atom.name.trim());
                    if rec_atom_type == "MMBBJ" {
                        model.membrane.push(atom_index as usize);
                    }

                    if active_restraints.contains(&res_id) {
                        match model.active_restraints.get_mut(&res_id) {
                            Some(atom_indexes) => {
                                atom_indexes.push(atom_index as usize);
                            },
                            None => {
                                model.active_restraints.insert(res_id.to_string(), vec![atom_index as usize]);
                            },
                        }
                    }

                    if passive_restraints.contains(&res_id) {
                        match model.passive_restraints.get_mut(&res_id) {
                            Some(atom_indexes) => {
                                atom_indexes.push(atom_index as usize);
                            },
                            None => {
                                model.passive_restraints.insert(res_id.to_string(), vec![atom_index as usize]);
                            },
                        }
                    }

                    let atom_name = atom.name.trim();
                    let mut atom_id = format!("{}-{}", residue.name.trim(), atom_name);

                    // Calculate AMBER type
                    let amber_type = match AMBER_TYPES.get(&*atom_id) {
                        Some(&amber) => amber,
                        _ => {
                            if atom_name == "H1" || atom_name == "H2" || atom_name == "H3" {
                                atom_id = format!("{}-H", residue.name.trim());
                                match AMBER_TYPES.get(&*atom_id) {
                                    Some(&amber) => amber,
                                    _ => panic!("DNA Error: Atom [{:?}] not supported", atom_id),
                                }
                            } else {
                                panic!("DNA Error: Atom [{:?}] not supported", atom_id);
                            }
                        },
                    };

                    // Assign electrostatics charge
                    let ele_charge = match ELE_CHARGES.get(&*atom_id) {
                        Some(&charge) => charge,
                        _ => {
                            match NT_ELE_CHARGES.get(&*atom_id) {
                                Some(&charge) => charge,
                                _ => panic!("DNA Error: Atom [{:?}] electrostatics charge not found", atom_id),
                            }
                        }
                    };
                    model.ele_charges.push(ele_charge);

                    // Assign VDW charge and radius
                    let vdw_charge = match VDW_CHARGES.get(&*amber_type) {
                        Some(&charge) => charge,
                        _ => panic!("DNA Error: Atom [{:?}] VDW charge not found", atom_id),
                    };
                    model.vdw_charges.push(vdw_charge);
                    let vdw_radius = match VDW_RADII.get(&*amber_type) {
                        Some(&radius) => radius,
                        _ => panic!("DNA Error: Atom [{:?}] VDW radius not found", atom_id),
                    };
                    model.vdw_radii.push(vdw_radius);

                    model.coordinates.push([atom.coord[0] as f64, atom.coord[1] as f64, atom.coord[2] as f64]);
                    atom_index += 1;
                }
            }
        }
        model
    }
}

pub struct DNA {
	pub potential: Vec<f64>,
    pub receptor: DNADockingModel,
    pub ligand: DNADockingModel,
    pub use_anm: bool,
}


impl<'a> DNA {

    pub fn new(receptor: Structure, rec_active_restraints: Vec<String>, rec_passive_restraints: Vec<String>,
            rec_nmodes: Vec<f64>, rec_num_anm: usize,
            ligand: Structure, lig_active_restraints: Vec<String>, lig_passive_restraints: Vec<String>,
            lig_nmodes: Vec<f64>, lig_num_anm: usize, use_anm: bool) -> Box<dyn Score + 'a> {
        let d = DNA {
            potential: Vec::with_capacity(168 * 168 * 20),
            receptor: DNADockingModel::new(&receptor, &rec_active_restraints, &rec_passive_restraints, &rec_nmodes, rec_num_anm),
            ligand: DNADockingModel::new(&ligand, &lig_active_restraints, &lig_passive_restraints, &lig_nmodes, lig_num_anm),
            use_anm,
        };
        Box::new(d)
    }
}

impl<'a> Score for DNA {

    fn energy(&self, translation: &[f64], rotation: &Quaternion,
        rec_nmodes: &[f64], lig_nmodes: &[f64]) -> f64 {

        // Clone receptor coordinates
        let mut receptor_coordinates: Vec<[f64; 3]> = self.receptor.coordinates.clone();
        let rec_num_atoms = receptor_coordinates.len();
        // Clone ligand coordinates
        let mut ligand_coordinates: Vec<[f64; 3]> = self.ligand.coordinates.clone();
        let lig_num_atoms = ligand_coordinates.len();

        // Get the proper ligand pose
        for (i_atom, coordinate) in ligand_coordinates.iter_mut().enumerate() {
            // First rotate
            let rotated_coordinate = rotation.rotate(coordinate.to_vec());
            // Then tranlate
            coordinate[0] = rotated_coordinate[0] + translation[0];
            coordinate[1] = rotated_coordinate[1] + translation[1];
            coordinate[2] = rotated_coordinate[2] + translation[2];
            // ANM
            if self.use_anm && self.ligand.num_anm > 0 {
                for i_nm in 0usize..self.ligand.num_anm {
                    // (num_anm, num_atoms, 3) -> 1d
                    // Endianness: i = i_nm * num_atoms * 3 + i_atom * 3 + coord
                    coordinate[0] += self.ligand.nmodes[i_nm * lig_num_atoms * 3 + i_atom * 3] * lig_nmodes[i_nm];
                    coordinate[1] += self.ligand.nmodes[i_nm * lig_num_atoms * 3 + i_atom * 3 + 1] * lig_nmodes[i_nm];
                    coordinate[2] += self.ligand.nmodes[i_nm * lig_num_atoms * 3 + i_atom * 3 + 2] * lig_nmodes[i_nm];
                }
            }
        }
        // Receptor only needs to use ANM
        for (i_atom, coordinate) in receptor_coordinates.iter_mut().enumerate() {
            // ANM
            if self.use_anm && self.receptor.num_anm > 0 {
                for i_nm in 0usize..self.receptor.num_anm {
                    // (num_anm, num_atoms, 3) -> 1d
                    // Endianness: i = i_nm * num_atoms * 3 + i_atom * 3 + coord
                    coordinate[0] += self.receptor.nmodes[i_nm * rec_num_atoms * 3 + i_atom * 3] * rec_nmodes[i_nm];
                    coordinate[1] += self.receptor.nmodes[i_nm * rec_num_atoms * 3 + i_atom * 3 + 1] * rec_nmodes[i_nm];
                    coordinate[2] += self.receptor.nmodes[i_nm * rec_num_atoms * 3 + i_atom * 3 + 2] * rec_nmodes[i_nm];
                }
            }
        }
        // Calculate scoring and interface
        let mut interface_receptor: Vec<usize> = vec![0; receptor_coordinates.len()];
        let mut interface_ligand: Vec<usize> = vec![0; ligand_coordinates.len()];

        let mut total_elec = 0.0;
        let mut total_vdw = 0.0;
        for (i, ra) in receptor_coordinates.iter().enumerate() {
            let x1 = ra[0];
            let y1 = ra[1];
            let z1 = ra[2];
            for (j, la) in ligand_coordinates.iter().enumerate() {
                let distance2 = (x1-la[0])*(x1-la[0]) + (y1-la[1])*(y1-la[1]) + (z1-la[2])*(z1-la[2]);

                // Electrostatics energy
                if distance2 <= ELEC_DIST_CUTOFF2 {
                    let mut atom_elec = self.receptor.ele_charges[i] * self.ligand.ele_charges[j] / distance2;
                    if atom_elec > ELEC_MAX_CUTOFF {
                        atom_elec = ELEC_MAX_CUTOFF;
                    }
                    if atom_elec < ELEC_MIN_CUTOFF {
                        atom_elec = ELEC_MIN_CUTOFF;
                    }
                    total_elec += atom_elec;
                }

                // Van der Waals energy
                if distance2 <= VDW_DIST_CUTOFF2 {
                    let vdw_energy = (self.receptor.vdw_charges[i] * self.ligand.vdw_charges[j]).sqrt();
                    let vdw_radius = self.receptor.vdw_radii[i] + self.ligand.vdw_radii[j];
                    let p6 = vdw_radius.powi(6) / distance2.powi(3);
                    let mut k = vdw_energy * (p6*p6 - 2.0 * p6);
                    if k > VDW_CUTOFF {
                        k = VDW_CUTOFF;
                    }
                    total_vdw += k;
                }

                // Interface calculation
                if distance2 <= INTERFACE_CUTOFF2 {
                    interface_receptor[i] = 1;
                    interface_ligand[j] = 1;
                }
            }
        }
        total_elec = total_elec * FACTOR / EPSILON;
        let score = (total_elec + total_vdw) * -1.0;

        // Bias the scoring depending on satisfied restraints
        let perc_receptor_restraints: f64 = satisfied_restraints(&interface_receptor,
            &self.receptor.active_restraints);
        let perc_ligand_restraints: f64 = satisfied_restraints(&interface_ligand,
            &self.ligand.active_restraints);
        // Take into account membrane intersection
        let mut membrane_penalty: f64 = 0.0;
        let intersection = membrane_intersection(&interface_receptor, &self.receptor.membrane);
        if intersection > 0.0 {
            membrane_penalty = MEMBRANE_PENALTY_SCORE * intersection;
        }

        score + perc_receptor_restraints * score + perc_ligand_restraints * score - membrane_penalty
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::env;
    use lib3dmol::parser;
    use crate::qt::Quaternion;

    #[test]
    fn test_1azp() {
        let cargo_path = match env::var("CARGO_MANIFEST_DIR") {
            Ok(val) => val,
            Err(_) => String::from("."),
        };
        let test_path: String = format!("{}/tests/1azp", cargo_path);

        let receptor_filename: String = format!("{}/1azp_receptor.pdb", test_path);
        let receptor = parser::read_pdb(&receptor_filename, "receptor");

        let ligand_filename: String = format!("{}/1azp_ligand.pdb", test_path);
        let ligand = parser::read_pdb(&ligand_filename, "ligand");

        let scoring = DNA::new(receptor, Vec::new(), Vec::new(), Vec::new(), 0, ligand, Vec::new(), Vec::new(), Vec::new(), 0, false);

        let translation = vec![0., 0., 0.];
        let rotation = Quaternion::default();
        let energy = scoring.energy(&translation, &rotation, &Vec::new(), &Vec::new());
        assert_eq!(energy, -364.88131078632557);
    }
}

