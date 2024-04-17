extern crate npyz;
extern crate serde;
extern crate serde_json;

use lightdock::constants::{
    DEFAULT_LIGHTDOCK_PREFIX, DEFAULT_LIG_NM_FILE, DEFAULT_REC_NM_FILE, DEFAULT_SEED,
};
use lightdock::dfire::DFIRE;
use lightdock::dna::DNA;
use lightdock::pydock::PYDOCK;
use lightdock::scoring::{Method, Score};
use lightdock::GSO;
use npyz::NpyFile;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use std::thread;

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
    let contents = fs::read_to_string(swarm_filename).expect("Error reading the input file");

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
    env_logger::init();
    // Parse command line
    let args: Vec<String> = env::args().collect();
    match args.len() {
        5 => {
            let setup_filename = &args[1];
            let swarm_filename = &args[2];
            let num_steps = &args[3];
            // parse the number
            let steps: u32 = match num_steps.parse() {
                Ok(n) => n,
                Err(_) => {
                    eprintln!("Error: steps argument must be a number");
                    return;
                }
            };
            let method_type = &args[4].to_lowercase();
            // parse the type
            let method = match &method_type[..] {
                "dfire" => Method::DFIRE,
                "dna" => Method::DNA,
                "pydock" => Method::PYDOCK,
                _ => {
                    eprintln!("Error: method not supported");
                    return;
                }
            };

            // Load setup
            let setup = read_setup_from_file(setup_filename).unwrap();

            // Simulation path
            let simulation_path = Path::new(setup_filename).parent().unwrap();

            simulate(
                simulation_path.to_str().unwrap(),
                &setup,
                swarm_filename,
                steps,
                method,
            );
        }
        _ => {
            println!(
                "Wrong command line. Usage: {} setup_filename swarm_filename steps method",
                args[0]
            );
        }
    }
}

fn parse_swarm_id(path: &Path) -> Option<i32> {
    path.file_name()
        .and_then(|s| s.to_str())
        .and_then(|s| s.strip_prefix("initial_positions_"))
        .and_then(|s| s.strip_suffix(".dat"))
        .and_then(|s| s.parse::<i32>().ok())
}

fn simulate(
    simulation_path: &str,
    setup: &SetupFile,
    swarm_filename: &str,
    steps: u32,
    method: Method,
) {
    let seed: u64 = match setup.seed {
        Some(seed) => seed,
        None => DEFAULT_SEED,
    };

    println!("Reading starting positions from {:?}", swarm_filename);
    let file_path = Path::new(swarm_filename);
    let swarm_id = parse_swarm_id(file_path).expect("Could not parse swarm from swarm filename");
    println!("Swarm ID {:?}", swarm_id);
    let swarm_directory = format!("swarm_{}", swarm_id);

    if !fs::metadata(&swarm_directory)
        .map(|m| m.is_dir())
        .unwrap_or(false)
    {
        panic!("Output directory does not exist for swarm {:?}", swarm_id);
    }

    println!("Writing to swarm dir {:?}", swarm_directory);
    let positions = parse_input_coordinates(swarm_filename);

    let receptor_filename = if simulation_path.is_empty() {
        format!("{}{}", DEFAULT_LIGHTDOCK_PREFIX, setup.receptor_pdb)
    } else {
        format!(
            "{}/{}{}",
            simulation_path, DEFAULT_LIGHTDOCK_PREFIX, setup.receptor_pdb
        )
    };
    // Parse receptor input PDB structure
    println!("Reading receptor input structure: {}", receptor_filename);
    let (receptor, _errors) =
        pdbtbx::open(&receptor_filename, pdbtbx::StrictnessLevel::Medium).unwrap();

    let ligand_filename = if simulation_path.is_empty() {
        format!("{}{}", DEFAULT_LIGHTDOCK_PREFIX, setup.ligand_pdb)
    } else {
        format!(
            "{}/{}{}",
            simulation_path, DEFAULT_LIGHTDOCK_PREFIX, setup.ligand_pdb
        )
    };
    // Parse ligand input PDB structure
    println!("Reading ligand input structure: {}", ligand_filename);
    let (ligand, _errors) =
        pdbtbx::open(&ligand_filename, pdbtbx::StrictnessLevel::Medium).unwrap();

    // Read ANM data if activated
    let mut rec_nm: Vec<f64> = Vec::new();
    let mut lig_nm: Vec<f64> = Vec::new();
    if setup.use_anm {
        if setup.anm_rec > 0 {
            let bytes = match std::fs::read(DEFAULT_REC_NM_FILE) {
                Ok(bytes) => bytes,
                Err(e) => {
                    panic!(
                        "Error reading receptor ANM file [{:?}]: {:?}",
                        DEFAULT_REC_NM_FILE,
                        e.to_string()
                    );
                }
            };
            let reader = NpyFile::new(&bytes[..]).unwrap();
            rec_nm = reader.into_vec::<f64>().unwrap();
            if rec_nm.len() != receptor.atom_count() * 3 * setup.anm_rec {
                panic!("Number of read ANM in receptor does not correspond to the number of atoms");
            }
        }
        if setup.anm_lig > 0 {
            let bytes = match std::fs::read(DEFAULT_LIG_NM_FILE) {
                Ok(bytes) => bytes,
                Err(e) => {
                    panic!(
                        "Error reading ligand ANM file [{:?}]: {:?}",
                        DEFAULT_LIG_NM_FILE,
                        e.to_string()
                    );
                }
            };
            let reader = NpyFile::new(&bytes[..]).unwrap();
            lig_nm = reader.into_vec::<f64>().unwrap();
            if lig_nm.len() != ligand.atom_count() * 3 * setup.anm_lig {
                panic!("Number of read ANM in ligand does not correspond to the number of atoms");
            }
        }
    }

    // Restraints
    let rec_active_restraints: Vec<String> = match &setup.receptor_restraints {
        Some(restraints) => restraints["active"].clone(),
        None => Vec::new(),
    };
    let rec_passive_restraints: Vec<String> = match &setup.receptor_restraints {
        Some(restraints) => restraints["passive"].clone(),
        None => Vec::new(),
    };
    let lig_active_restraints: Vec<String> = match &setup.ligand_restraints {
        Some(restraints) => restraints["active"].clone(),
        None => Vec::new(),
    };
    let lig_passive_restraints: Vec<String> = match &setup.ligand_restraints {
        Some(restraints) => restraints["passive"].clone(),
        None => Vec::new(),
    };

    // Scoring function
    println!("Loading {:?} scoring function", method);
    let scoring = match method {
        Method::DFIRE => DFIRE::new(
            receptor,
            rec_active_restraints,
            rec_passive_restraints,
            rec_nm,
            setup.anm_rec,
            ligand,
            lig_active_restraints,
            lig_passive_restraints,
            lig_nm,
            setup.anm_lig,
            setup.use_anm,
        ) as Box<dyn Score>,
        Method::DNA => DNA::new(
            receptor,
            rec_active_restraints,
            rec_passive_restraints,
            rec_nm,
            setup.anm_rec,
            ligand,
            lig_active_restraints,
            lig_passive_restraints,
            lig_nm,
            setup.anm_lig,
            setup.use_anm,
        ) as Box<dyn Score>,
        Method::PYDOCK => PYDOCK::new(
            receptor,
            rec_active_restraints,
            rec_passive_restraints,
            rec_nm,
            setup.anm_rec,
            ligand,
            lig_active_restraints,
            lig_passive_restraints,
            lig_nm,
            setup.anm_lig,
            setup.use_anm,
        ) as Box<dyn Score>,
    };

    // Glowworm Swarm Optimization algorithm
    println!("Creating GSO with {} glowworms", positions.len());
    let mut gso = GSO::new(
        &positions,
        seed,
        &scoring,
        setup.use_anm,
        setup.anm_rec,
        setup.anm_lig,
        swarm_directory,
    );

    // Simulate for the given steps
    println!("Starting optimization ({} steps)", steps);
    gso.run(steps);
}
