#[macro_use]
extern crate lazy_static;
extern crate rand;

pub mod constants;
pub mod dfire;
pub mod dna;
pub mod glowworm;
pub mod pydock;
pub mod qt;
pub mod scoring;
pub mod swarm;

use log::info;
use rand::rngs::StdRng;
use rand::SeedableRng;
use scoring::Score;
use swarm::Swarm;

pub struct GSO<'a> {
    pub swarm: Swarm<'a>,
    pub rng: StdRng,
    pub output_directory: String,
}

impl<'a> GSO<'a> {
    pub fn new(
        positions: &[Vec<f64>],
        seed: u64,
        scoring: &'a Box<dyn Score>,
        use_anm: bool,
        rec_num_anm: usize,
        lig_num_anm: usize,
        output_directory: String,
    ) -> Self {
        let mut gso = GSO {
            swarm: Swarm::new(),
            rng: SeedableRng::seed_from_u64(seed),
            output_directory,
        };
        gso.swarm
            .add_glowworms(positions, scoring, use_anm, rec_num_anm, lig_num_anm);
        gso
    }

    pub fn run(&mut self, steps: u32) {
        for step in 1..steps + 1 {
            info!("Step {}", step);
            self.swarm.update_luciferin();
            self.swarm.movement_phase(&mut self.rng);
            if step % 10 == 0 || step == 1 {
                match self.swarm.save(step, &self.output_directory) {
                    Ok(ok) => ok,
                    Err(why) => panic!("Error saving GSO output: {:?}", why),
                }
            }
        }
    }
}
