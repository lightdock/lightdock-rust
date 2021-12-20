use lib3dmol::structures::Structure;
use std::collections::HashMap;

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
        _ => { panic!("Residue name not supported in DFIRE scoring function") }
    }
}

pub struct DockingModel {
    pub atoms: Vec<usize>,
    pub coordinates: Vec<[f64; 3]>,
    pub membrane: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>,
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub nmodes: Vec<f64>,
    pub num_anm: usize,
}

#[derive(Debug)]
pub enum Method {
    DFIRE,
}

pub trait Score {
    fn energy(&self, receptor: &DockingModel, ligand: &DockingModel,
        receptor_coordinates: &[[f64; 3]], ligand_coordinates: &[[f64; 3]],
        interface_receptor: &mut Vec<usize>, interface_ligand: &mut Vec<usize>) -> f64;

    fn get_docking_model(&self, structure: &Structure,
        active_restraints: &[String], passive_restraints: &[String],
        nmodes: &Vec<f64>, num_anm: usize) -> DockingModel;
}
