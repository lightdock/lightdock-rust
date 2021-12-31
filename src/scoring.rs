use std::collections::HashMap;
use super::qt::Quaternion;


#[derive(Debug)]
pub enum Method {
    DFIRE,
    DNA,
}

pub trait Score {
    fn energy(&self, translation: &[f64], rotation: &Quaternion,
        rec_nmodes: &[f64], lig_nmodes: &[f64]) -> f64;
}

pub fn satisfied_restraints(interface: &[usize], restraints: &HashMap<String, Vec<usize>>) -> f64 {
    // Calculate the percentage of satisfied restraints
    if restraints.is_empty() {
        return 0.0
    }
    let mut num_residues = 0;
    for (_k, atom_indexes) in restraints.iter() {
        for &i in atom_indexes.iter() {
            if interface[i] == 1 {
                num_residues += 1;
                break;
            }
        }
    }
    num_residues as f64 / restraints.len() as f64
}

pub fn membrane_intersection(interface: &[usize], membrane: &[usize]) -> f64 {
    if membrane.is_empty() {
        return 0.0
    }
    let mut num_beads = 0;
    for &i_bead in membrane.iter() {
        num_beads += interface[i_bead];
    }
    num_beads as f64 / membrane.len() as f64
}