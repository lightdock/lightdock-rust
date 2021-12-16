use std::f64;
use std::collections::HashMap;
use super::scoring::{Score, DockingModel};
use super::qt::Quaternion;
use super::constants::{DEFAULT_TRANSLATION_STEP, DEFAULT_ROTATION_STEP, DEFAULT_NMODES_STEP, MEMBRANE_PENALTY_SCORE};


pub struct Glowworm<'a> {
    pub id: u32,
    pub translation: Vec<f64>,
    pub rotation: Quaternion,
    pub rec_nmodes: Vec<f64>,
    pub lig_nmodes: Vec<f64>,
    pub scoring_function: &'a dyn Score,
    pub receptor: &'a DockingModel,
    pub ligand: &'a DockingModel,
    pub rho: f64,
    pub gamma: f64,
    pub beta: f64,
    pub luciferin: f64,
    pub vision_range: f64,
    pub max_vision_range: f64,
    pub max_neighbors: u32,
    pub neighbors: Vec<u32>,
    pub probabilities: Vec<f64>,
    pub scoring: f64,
    pub moved: bool,
    pub step: u32,
    pub use_anm: bool,
}

impl<'a> Glowworm<'a> {
     pub fn new(id: u32, translation:Vec<f64>, rotation:Quaternion, 
        rec_nmodes: Vec<f64>, lig_nmodes: Vec<f64>, scoring_function: &'a dyn Score,
        receptor: &'a DockingModel, ligand: &'a DockingModel, use_anm: bool) -> Self {
        Glowworm {
            id,
            translation,
            rotation,
            rec_nmodes,
            lig_nmodes,
            scoring_function,
            receptor,
            ligand,
            rho: 0.5,
            gamma: 0.4,
            beta: 0.08,
            luciferin: 5.0,
            vision_range: 0.2,
            max_vision_range: 5.0,
            max_neighbors: 5,
            neighbors: Vec::new(),
            probabilities: Vec::new(),
            scoring: 0.0,
            moved: false,
            step: 0,
            use_anm,
        }
    }

    pub fn compute_luciferin(&mut self) {
        if self.moved || self.step == 0 {
            let mut receptor_coordinates: Vec<[f64; 3]> = self.receptor.coordinates.clone();
            let rec_num_atoms = receptor_coordinates.len();
            let mut ligand_coordinates: Vec<[f64; 3]> = self.ligand.coordinates.clone();
            let lig_num_atoms = ligand_coordinates.len();
            // Get the proper ligand pose
            for (i_atom, coordinate) in ligand_coordinates.iter_mut().enumerate() {
                // First rotate
                let rotated_coordinate = self.rotation.rotate(coordinate.to_vec());
                // Then tranlate
                coordinate[0] = rotated_coordinate[0] + self.translation[0];
                coordinate[1] = rotated_coordinate[1] + self.translation[1];
                coordinate[2] = rotated_coordinate[2] + self.translation[2];
                // ANM
                if self.use_anm && self.ligand.num_anm > 0 {
                    for i_nm in 0usize..self.ligand.num_anm {
                        // (num_anm, num_atoms, 3) -> 1d
                        // Endianness: i = i_nm * num_atoms * 3 + i_atom * 3 + coord
                        coordinate[0] += self.ligand.nmodes[i_nm * lig_num_atoms * 3 + i_atom * 3] * self.lig_nmodes[i_nm];
                        coordinate[1] += self.ligand.nmodes[i_nm * lig_num_atoms * 3 + i_atom * 3 + 1] * self.lig_nmodes[i_nm];
                        coordinate[2] += self.ligand.nmodes[i_nm * lig_num_atoms * 3 + i_atom * 3 + 2] * self.lig_nmodes[i_nm];
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
                        coordinate[0] += self.receptor.nmodes[i_nm * rec_num_atoms * 3 + i_atom * 3] * self.rec_nmodes[i_nm];
                        coordinate[1] += self.receptor.nmodes[i_nm * rec_num_atoms * 3 + i_atom * 3 + 1] * self.rec_nmodes[i_nm];
                        coordinate[2] += self.receptor.nmodes[i_nm * rec_num_atoms * 3 + i_atom * 3 + 2] * self.rec_nmodes[i_nm];
                    }
                }
            }
            // Calculate scoring and interface
            let mut interface_receptor: Vec<usize> = vec![0; receptor_coordinates.len()];
            let mut interface_ligand: Vec<usize> = vec![0; ligand_coordinates.len()];
            let energy = self.scoring_function.energy(self.receptor, self.ligand, 
                &receptor_coordinates, &ligand_coordinates, 
                &mut interface_receptor, &mut interface_ligand);
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

            self.scoring = energy + perc_receptor_restraints * energy 
                + perc_ligand_restraints * energy - membrane_penalty;
        }
        self.luciferin = (1.0 - self.rho) * self.luciferin + self.gamma * self.scoring;
        self.step += 1;
    }

    pub fn distance(&mut self, other: &Glowworm) -> f64 {
        let x1 = self.translation[0];
        let x2 = other.translation[0];
        let y1 = self.translation[1];
        let y2 = other.translation[1];
        let z1 = self.translation[2];
        let z2 = other.translation[2];
        ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)).sqrt()
    }

    pub fn is_neighbor(&mut self, other: &Glowworm) -> bool {
        if self.id != other.id && self.luciferin < other.luciferin {
            return self.distance(other) < self.vision_range;
        }
        false
    }

    pub fn update_vision_range(&mut self) {
        self.vision_range = (self.max_vision_range).min(
            (0_f64).max(self.vision_range + self.beta * f64::from(self.max_neighbors as i32 - 
                        (self.neighbors.len() as i32))));
    }

    pub fn compute_probability_moving_toward_neighbor(&mut self, luciferins: &[f64]) {
        self.probabilities = Vec::new();

        let mut total_sum: f64 = 0.0;
        let mut difference: f64;
        for neighbor_id in &self.neighbors {
            difference = luciferins[*neighbor_id as usize] - self.luciferin;
            self.probabilities.push(difference);
            total_sum += difference;
        }

        for i in 0..self.neighbors.len(){
            self.probabilities[i] /= total_sum;
        }
    }


    pub fn select_random_neighbor(&mut self, random_number:f64) -> u32 {
        if self.neighbors.is_empty() {
            return self.id;
        }

        let mut sum_probabilities:f64 = 0.0;
        let mut i:usize = 0;
        while sum_probabilities < random_number {
            sum_probabilities += self.probabilities[i];
            i += 1;
        }
        self.neighbors[i-1]
    }

    pub fn move_towards(&mut self, other_id:u32, other_position:&[f64], 
        other_rotation:&Quaternion, other_anm_rec:&[f64], other_anm_lig:&[f64]) {
        self.moved = self.id != other_id;
        if self.id != other_id {
            // Translation component
            let mut delta_x:Vec<f64> = vec![other_position[0] - self.translation[0],
                                            other_position[1] - self.translation[1],
                                            other_position[2] - self.translation[2]];
            let norm:f64 = (delta_x[0]*delta_x[0] + delta_x[1]*delta_x[1] + delta_x[2]*delta_x[2]).sqrt();
            let coef:f64 = DEFAULT_TRANSLATION_STEP / norm;
            delta_x[0] *= coef;
            delta_x[1] *= coef;
            delta_x[2] *= coef;
            self.translation[0] += delta_x[0];
            self.translation[1] += delta_x[1];
            self.translation[2] += delta_x[2];
            
            // Rotation component
            self.rotation = self.rotation.slerp(other_rotation, DEFAULT_ROTATION_STEP);
            
            // ANM component
            if self.use_anm && self.receptor.num_anm > 0 {
                let mut delta_anm:Vec<f64> = Vec::new();
                let mut cum_norm: f64 = 0.0;
                for i in 0..self.rec_nmodes.len() {
                    let diff = other_anm_rec[i] - self.rec_nmodes[i];
                    delta_anm.push(diff);
                    cum_norm += diff*diff
                }
                let anm_rec_norm:f64 = cum_norm.sqrt();
                let anm_rec_coef:f64 = DEFAULT_NMODES_STEP / anm_rec_norm;
                for i in 0..self.rec_nmodes.len() {
                    delta_anm[i] *= anm_rec_coef;
                    self.rec_nmodes[i] += delta_anm[i];
                }
            }
            if self.use_anm && self.ligand.num_anm > 0 {
                let mut delta_anm:Vec<f64> = Vec::new();
                let mut cum_norm: f64 = 0.0;
                for i in 0..self.lig_nmodes.len() {
                    let diff = other_anm_lig[i] - self.lig_nmodes[i];
                    delta_anm.push(diff);
                    cum_norm += diff*diff
                }
                let anm_lig_norm:f64 = cum_norm.sqrt();
                let anm_lig_coef:f64 = DEFAULT_NMODES_STEP / anm_lig_norm;
                for i in 0..self.lig_nmodes.len() {
                    delta_anm[i] *= anm_lig_coef;
                    self.lig_nmodes[i] += delta_anm[i];
                }
            }
        }
    }
}

pub fn distance(one: &Glowworm, two: &Glowworm) -> f64 {
    // Calculate the distance between two glowworms using their translation vector
    let x1 = one.translation[0];
    let x2 = two.translation[0];
    let y1 = one.translation[1];
    let y2 = two.translation[1];
    let z1 = one.translation[2];
    let z2 = two.translation[2];
    ((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)).sqrt()
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
