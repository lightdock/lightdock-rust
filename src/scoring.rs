use lib3dmol::structures::Structure;
use std::collections::HashMap;


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
    DNA,
}

pub trait Score {
    fn energy(&self, receptor: &DockingModel, ligand: &DockingModel,
        receptor_coordinates: &[[f64; 3]], ligand_coordinates: &[[f64; 3]],
        interface_receptor: &mut Vec<usize>, interface_ligand: &mut Vec<usize>) -> f64;

    fn get_docking_model(&self, structure: &Structure,
        active_restraints: &[String], passive_restraints: &[String],
        nmodes: &Vec<f64>, num_anm: usize) -> DockingModel;
}
