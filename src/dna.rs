use lib3dmol::structures::Structure;
use std::collections::HashMap;
use super::qt::Quaternion;
use super::constants::{INTERFACE_CUTOFF, MEMBRANE_PENALTY_SCORE};
use super::scoring::{Score, satisfied_restraints, membrane_intersection};

macro_rules! hashmap {
    ($( $key: expr => $val: expr ),*) => {{
         let mut map = ::std::collections::HashMap::new();
         $( map.insert($key, $val); )*
         map
    }}
}

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
}


pub struct DNADockingModel {
    pub atoms: Vec<usize>,
    pub coordinates: Vec<[f64; 3]>,
    pub membrane: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub num_anm: usize,
    pub nmodes: Vec<f64>,
    pub amber_types: Vec<String>,
}

impl<'a> DNADockingModel {

    fn new(structure: &'a Structure, active_restraints: &'a [String], passive_restraints: &'a [String],
        nmodes: &'a Vec<f64>, num_anm: usize) -> DNADockingModel {
        let mut model = DNADockingModel {
            atoms: Vec::new(),
            coordinates: Vec::new(),
            membrane: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
            nmodes: nmodes.clone(),
            num_anm: num_anm,
            amber_types: Vec::new(),
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
                    model.amber_types.push(amber_type.to_string());
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

    fn energy(&self, translation: &Vec<f64>, rotation: &Quaternion,
        rec_nmodes: &Vec<f64>, lig_nmodes: &Vec<f64>) -> f64 {
        let mut score: f64 = 0.0;

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

        for (i, ra) in receptor_coordinates.iter().enumerate() {
            let x1 = ra[0];
            let y1 = ra[1];
            let z1 = ra[2];
            for (j, la) in ligand_coordinates.iter().enumerate() {
                let dist = (x1-la[0])*(x1-la[0]) + (y1-la[1])*(y1-la[1]) + (z1-la[2])*(z1-la[2]);
                if dist <= 225. {
                    let d = dist.sqrt()*2.0 - 1.0;
                    score += 1.0;
                    if d <= INTERFACE_CUTOFF {
                        interface_receptor[i] = 1;
                        interface_ligand[j] = 1;
                    }
                }
            }
        }

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

    #[test]
    fn test_read_potentials() {
        let mut scoring = DNA {
            potential: Vec::with_capacity(168 * 168 * 20),
        };
        scoring.load_potentials();
        assert_eq!(scoring.potential[0], 10.0);
        assert_eq!(scoring.potential[2], -0.624030868);
        assert_eq!(scoring.potential[4998], -0.0458685914);
        assert_eq!(scoring.potential[168*168*20-1], 0.0);
    }

    #[test]
    fn test_2oob() {
        let cargo_path = match env::var("CARGO_MANIFEST_DIR") {
            Ok(val) => val,
            Err(_) => String::from("."),
        };
        let test_path: String = format!("{}/tests/2oob", cargo_path);

        let scoring = DNA::new();

        let receptor_filename: String = format!("{}/2oob_receptor.pdb", test_path);
        let receptor = parser::read_pdb(&receptor_filename, "receptor");
        let receptor_model = scoring.get_docking_model(&receptor, &Vec::new(), &Vec::new(), &Vec::new(), 0);

        let ligand_filename: String = format!("{}/2oob_ligand.pdb", test_path);
        let ligand = parser::read_pdb(&ligand_filename, "ligand");
        let ligand_model = scoring.get_docking_model(&ligand, &Vec::new(), &Vec::new(), &Vec::new(), 0);

        let mut receptor_coordinates: Vec<[f64; 3]> = receptor_model.coordinates.clone();
        let mut ligand_coordinates: Vec<[f64; 3]> = ligand_model.coordinates.clone();
        let mut interface_receptor: Vec<usize> = vec![0; receptor_coordinates.len()];
        let mut interface_ligand: Vec<usize> = vec![0; ligand_coordinates.len()];
        let energy = scoring.energy(&receptor_model, &ligand_model, &receptor_coordinates, &ligand_coordinates,
                                    &mut interface_receptor, &mut interface_ligand);
        assert_eq!(energy, 16.7540569503498);
    }
}
