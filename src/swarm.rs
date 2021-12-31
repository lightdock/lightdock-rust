use super::glowworm::Glowworm;
use super::glowworm::distance;
use super::scoring::Score;
use super::qt::Quaternion;
use rand::Rng;
use std::fs::File;
use std::io::{Write, Error};


pub struct Swarm<'a> {
	pub glowworms: Vec<Glowworm<'a>>,
}

impl<'a> Default for Swarm<'a> {
    fn default() -> Self {
        Swarm::new()
    }
}

impl<'a> Swarm<'a> {
    pub fn new() -> Self {
    	Swarm {
    		glowworms: Vec::new(),
    	}
    }

    pub fn add_glowworms(&mut self, positions: &[Vec<f64>], scoring: &'a Box<dyn Score>, use_anm: bool,
        rec_num_anm: usize, lig_num_anm: usize) {
    	for (i, position) in positions.iter().enumerate() {
            // Translation component
            let translation = vec![position[0], position[1], position[2]];
            // Rotation component
            let rotation = Quaternion::new(position[3], position[4], position[5], position[6]);
            // ANM for receptor
            let mut rec_nmodes: Vec<f64> = Vec::new();
            if use_anm && rec_num_anm > 0 {
                for j in 7..7+rec_num_anm {
                    rec_nmodes.push(positions[i][j as usize]);
                }
            }
            // ANM for ligand
            let mut lig_nmodes: Vec<f64> = Vec::new();
            if use_anm && lig_num_anm > 0 {
                for j in 7+rec_num_anm..positions[i].len() {
                    lig_nmodes.push(positions[i][j as usize]);
                }
            }
            let glowworm = Glowworm::new(i as u32, translation, rotation, rec_nmodes, lig_nmodes, scoring, use_anm);
    		self.glowworms.push(glowworm);
    	}
    }

    pub fn update_luciferin(&mut self) {
    	for glowworm in self.glowworms.iter_mut() {
    		glowworm.compute_luciferin();
    	}
    }

    pub fn movement_phase(&mut self, rng: &mut rand::prelude::StdRng) {
    	// Save original positions
    	let mut positions: Vec<Vec<f64>> = Vec::new();
        let mut rotations: Vec<Quaternion> = Vec::new();
        let mut anm_recs: Vec<Vec<f64>> = Vec::new();
        let mut anm_ligs: Vec<Vec<f64>> = Vec::new();
    	for glowworm in self.glowworms.iter(){
    		positions.push(glowworm.translation.clone());
            rotations.push(glowworm.rotation);
            anm_recs.push(glowworm.rec_nmodes.clone());
            anm_ligs.push(glowworm.lig_nmodes.clone());
    	}

    	// First search for each glowworm's neighbors
    	let mut neighbors: Vec<Vec<u32>> = Vec::new();
    	for i in 0..self.glowworms.len() {
    		let mut this_neighbors = Vec::new();
    		let g1 = &self.glowworms[i];
        	for j in 0..self.glowworms.len() {
        		if i != j {
        			let g2 = &self.glowworms[j];
        			if g1.luciferin < g2.luciferin {
        				let distance = distance(g1, g2);
            			if distance < g1.vision_range {
            				this_neighbors.push(g2.id);
            			}   
            		}
        		}
        	}
        	neighbors.push(this_neighbors);
    	}
    	
    	// Second compute probability moving towards the neighbor
    	let mut luciferins = Vec::new();
    	for glowworm in self.glowworms.iter_mut() {
    		luciferins.push(glowworm.luciferin);
    	}
    	for i in 0..self.glowworms.len() {
    		let glowworm = &mut self.glowworms[i];
    		glowworm.neighbors = neighbors[i].clone();
    		glowworm.compute_probability_moving_toward_neighbor(&luciferins);
    	}
    	
    	// Finally move to the selected position
        for i in 0..self.glowworms.len() {
        	let glowworm = &mut self.glowworms[i];
            let neighbor_id = glowworm.select_random_neighbor(rng.gen::<f64>());
            let position = &positions[neighbor_id as usize];
            let rotation = &rotations[neighbor_id as usize];
            let anm_rec = &anm_recs[neighbor_id as usize];
            let anm_lig = &anm_ligs[neighbor_id as usize];
            glowworm.move_towards(neighbor_id, position, rotation, anm_rec, anm_lig);
            glowworm.update_vision_range();
        }
    }

    pub fn save(&mut self, step: u32) -> Result<(), Error> {
        let path = format!("gso_{:?}.out", step);
        let mut output = File::create(path)?;
        writeln!(output, "#Coordinates  RecID  LigID  Luciferin  Neighbor's number  Vision Range  Scoring")?;
        for glowworm in self.glowworms.iter() {
            write!(output, "({:.7}, {:.7}, {:.7}, {:.7}, {:.7}, {:.7}, {:.7}", 
                glowworm.translation[0], glowworm.translation[1], glowworm.translation[2],
                glowworm.rotation.w, glowworm.rotation.x, glowworm.rotation.y, glowworm.rotation.z)?;
            if glowworm.use_anm && !glowworm.rec_nmodes.is_empty() {
                for i in 0..glowworm.rec_nmodes.len() {
                    write!(output, ", {:.7}", glowworm.rec_nmodes[i])?;
                }
            }
            if glowworm.use_anm && !glowworm.lig_nmodes.is_empty() {
                for i in 0..glowworm.lig_nmodes.len() {
                    write!(output, ", {:.7}", glowworm.lig_nmodes[i])?;
                }
            }
            writeln!(output, ")    0    0   {:.8}  {:?} {:.3} {:.8}",
                glowworm.luciferin, glowworm.neighbors.len(), glowworm.vision_range, glowworm.scoring)?;
        }
        Ok(())
    }
}
