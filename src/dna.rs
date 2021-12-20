use lib3dmol::structures::Structure;
use std::collections::HashMap;

use super::scoring::{DockingModel, Score};

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

pub struct DNA {
}

impl Score for DNA {

    fn get_docking_model(&self, structure: &Structure,
        active_restraints: &[String], passive_restraints: &[String],
        nmodes: &Vec<f64>, num_anm: usize) -> DockingModel {
        let mut model = DockingModel {
            atoms: Vec::new(),
            coordinates: Vec::new(),
            membrane: Vec::new(),
            active_restraints: HashMap::new(),
            passive_restraints: HashMap::new(),
            nmodes: nmodes.clone(),
            num_anm: num_anm,
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

                    model.coordinates.push([atom.coord[0] as f64, atom.coord[1] as f64, atom.coord[2] as f64]);
                    atom_index += 1;
                }
            }
        }
        model
    }

    fn energy(&self, receptor: &DockingModel, ligand: &DockingModel,
        receptor_coordinates: &[[f64; 3]], ligand_coordinates: &[[f64; 3]],
        interface_receptor: &mut Vec<usize>, interface_ligand: &mut Vec<usize>) -> f64 {
        0.0
    }
}

impl DNA {
    pub fn new() -> Box<dyn Score> {
        Box::new(DNA{})
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_potentials() {
        let scoring = DNA::new();
    }
}
