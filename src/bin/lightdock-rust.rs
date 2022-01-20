extern crate serde;
extern crate serde_json;
extern crate npy;

use lightdock::GSO;
use lightdock::constants::{DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_SEED, DEFAULT_REC_NM_FILE, DEFAULT_LIG_NM_FILE};
use lightdock::scoring::{Score, Method};
use lightdock::dfire::DFIRE;
use lightdock::dna::DNA;
use std::env;
use std::fs;
use serde::{Serialize, Deserialize};
use std::error::Error;
use std::fs::File;
use std::io::{Read, BufReader};
use std::path::Path;
use std::collections::HashMap;
use std::thread;
use lib3dmol::parser;
use npy::NpyData;

// Use 8MB as binary stack
const STACK_SIZE: usize = 8 * 1024 * 1024;

#[derive(Serialize, Deserialize, Debug)]
struct SetupFile {
    seed: Option<u64>, 
    anm_seed: u64, 
    ftdock_file: Option<String>, 
    noh: bool, 
    anm_rec: usize, 
    anm_lig: usize, 
    swarms: u32, 
    starting_points_seed: u32, 
    verbose_parser: bool, 
    noxt: bool, 
    now: bool,
    restraints: Option<String>, 
    use_anm: bool, 
    glowworms: u32, 
    membrane: bool, 
    receptor_pdb: String, 
    ligand_pdb: String,
    receptor_restraints: Option<HashMap<String, Vec<String>>>,
    ligand_restraints: Option<HashMap<String, Vec<String>>>,
}

fn read_setup_from_file<P: AsRef<Path>>(path: P) -> Result<SetupFile, Box<dyn Error>> {
    // Open the file in read-only mode with buffer.
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    // Read the JSON contents of the file as an instance of `SetupFile`.
    let u = serde_json::from_reader(reader)?;
    // Return the `SetupFile`.
    Ok(u)
}

fn parse_input_coordinates(swarm_filename: &str) -> Vec<Vec<f64>> {
    // Parse swarm filename content
    let contents = fs::read_to_string(swarm_filename)
        .expect("Error reading the input file");

    let mut positions: Vec<Vec<f64>> = Vec::new();
    for s in contents.lines() {
        let vector_raw: String = String::from(s);
        let vector: Vec<&str> = vector_raw.split(' ').collect();
        let mut position: Vec<f64> = Vec::new();
        for pos in vector.iter() {         
            position.push(pos.trim().parse::<f64>().unwrap());
        }
        positions.push(position);
    }
    positions
}

fn main() {
    // Spawn thread with explicit stack size
    let child = thread::Builder::new()
        .stack_size(STACK_SIZE)
        .spawn(run)
        .unwrap();

    // Wait for thread to join
    child.join().unwrap();
}

fn run() {
    // Parse command line
    let args: Vec<String> = env::args().collect();
    match args.len() {
        5 => {
            let setup_filename = &args[1];
            let swarm_filename = &args[2];
            let num_steps = &args[3];
            // parse the number
            let steps: u32 = match num_steps.parse() {
                Ok(n) => {
                    n
                },
                Err(_) => {
                    eprintln!("Error: steps argument must be a number");
                    return;
                },
            };
            let method_type = &args[4].to_lowercase();
            // parse the type
            let method = match &method_type[..] {
                "dfire" => Method::DFIRE,
                "dna" => Method::DNA,
                _ => {
                    eprintln!("Error: method not supported");
                    return;
                },
            };

            // Load setup
            let setup = read_setup_from_file(setup_filename).unwrap();

            simulate(&setup, swarm_filename, steps, method);
        }
        _ => {
            println!("Wrong command line. Usage: {} setup_filename swarm_filename steps method", args[0]);
        }
    }
}

fn simulate(setup: &SetupFile, swarm_filename: &str, steps: u32, method: Method) {

    let seed:u64 = match setup.seed {
        Some(seed) => {
            seed
        },
        None => {
            DEFAULT_SEED
        },
    };

    println!("Reading starting positions from {:?}", swarm_filename);
    let positions = parse_input_coordinates(swarm_filename);

    // Parse receptor input PDB structure
    let receptor_filename: String = format!("{}{}", DEFAULT_LIGHTDOCK_PREFIX, setup.receptor_pdb);
    println!("Reading receptor input structure: {}", receptor_filename);
    let receptor = parser::read_pdb(&receptor_filename, "receptor");

    // Parse ligand input PDB structure
    let ligand_filename: String = format!("{}{}", DEFAULT_LIGHTDOCK_PREFIX, setup.ligand_pdb);
    println!("Reading ligand input structure: {}", ligand_filename);
    let ligand = parser::read_pdb(&ligand_filename, "ligand");

    // Read ANM data if activated
    let mut rec_nm: Vec<f64> = Vec::new();
    let mut lig_nm: Vec<f64> = Vec::new();
    if setup.use_anm {
        let mut buf = vec![];
        if setup.anm_rec > 0 {
            std::fs::File::open(DEFAULT_REC_NM_FILE).unwrap().read_to_end(&mut buf).unwrap();
            rec_nm = NpyData::from_bytes(&buf).unwrap().to_vec();
            if rec_nm.len() != receptor.get_atom_number() as usize * 3 * setup.anm_rec {
                panic!("Number of read ANM in receptor does not correspond to the number of atoms");
            }
        }
        if setup.anm_lig > 0 {
            buf = vec![];
            std::fs::File::open(DEFAULT_LIG_NM_FILE).unwrap().read_to_end(&mut buf).unwrap();
            lig_nm = NpyData::from_bytes(&buf).unwrap().to_vec();
            if lig_nm.len() != ligand.get_atom_number() as usize * 3 * setup.anm_lig {
                panic!("Number of read ANM in ligand does not correspond to the number of atoms");
            }
        }
    }

    // Restraints
    let rec_active_restraints: Vec<String> = match &setup.receptor_restraints {
        Some(restraints) => { restraints["active"].clone() },
        None => { Vec::new() },
    };
    let rec_passive_restraints: Vec<String> = match &setup.receptor_restraints {
        Some(restraints) => { restraints["passive"].clone() },
        None => { Vec::new() },
    };
    let lig_active_restraints: Vec<String> = match &setup.ligand_restraints {
        Some(restraints) => { restraints["active"].clone() },
        None => { Vec::new() },
    };
    let lig_passive_restraints: Vec<String> = match &setup.ligand_restraints {
        Some(restraints) => { restraints["passive"].clone() },
        None => { Vec::new() },
    };

    // Scoring function
    println!("Loading {:?} scoring function", method);
    let scoring = match method {
        Method::DFIRE => DFIRE::new(receptor, rec_active_restraints, rec_passive_restraints, rec_nm, setup.anm_rec,
            ligand, lig_active_restraints, lig_passive_restraints, lig_nm, setup.anm_lig, setup.use_anm) as Box<dyn Score>,
        Method::DNA => DNA::new(receptor, rec_active_restraints, rec_passive_restraints, rec_nm, setup.anm_rec,
            ligand, lig_active_restraints, lig_passive_restraints, lig_nm, setup.anm_lig, setup.use_anm) as Box<dyn Score>,
    };

    // Glowworm Swarm Optimization algorithm
    println!("Creating GSO with {} glowworms", positions.len());
    let mut gso = GSO::new(&positions, seed, &scoring, setup.use_anm, setup.anm_rec, setup.anm_lig);

    // Simulate for the given steps
    println!("Starting optimization ({} steps)", steps);
    gso.run(steps);
}
