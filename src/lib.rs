#[macro_use]

extern crate lazy_static;
extern crate rand;

pub mod glowworm;
pub mod swarm;
pub mod qt;
pub mod constants;
pub mod scoring;
pub mod dfire;
pub mod dna;

use swarm::Swarm;
use rand::SeedableRng;
use rand::rngs::StdRng;
use scoring::Score;


pub struct GSO<'a> {
    pub swarm: Swarm<'a>,
    pub rng: StdRng,
}


impl<'a> GSO<'a> {
    pub fn new(positions: &[Vec<f64>], seed: u64, scoring: &'a Box<dyn Score>, use_anm: bool,
        rec_num_anm: usize, lig_num_anm: usize) -> Self {
        let mut gso = GSO {
            swarm: Swarm::new(),
            rng: SeedableRng::seed_from_u64(seed),
        };
        gso.swarm.add_glowworms(positions, scoring, use_anm, rec_num_anm, lig_num_anm);
        gso
    }

    pub fn run(&mut self, steps: u32) {
        for step in 1..steps+1 {
            println!("Step {}", step);
            self.swarm.update_luciferin();
            self.swarm.movement_phase(&mut self.rng);
            if step % 10 == 0 || step == 1 {
                match self.swarm.save(step) {
                    Ok(ok) => ok,
                    Err(why) => panic!("Error saving GSO output: {:?}", why),
                }
            }
        }
    }
}
