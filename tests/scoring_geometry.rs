//! Tests for the ligand pose and the interface calculated by the scoring functions
use lightdock::dfire::DFIRE;
use lightdock::dna::DNA;
use lightdock::pydock::PYDOCK;
use lightdock::qt::Quaternion;
use std::fs;

fn open(path: &str) -> pdbtbx::PDB {
    pdbtbx::open(path, pdbtbx::StrictnessLevel::Medium)
        .unwrap()
        .0
}

/// Writes a structure made of a single alpha carbon at the origin
fn single_atom(name: &str, chain: char) -> String {
    let path = format!("{}/{}.pdb", env!("CARGO_TARGET_TMPDIR"), name);
    let line = format!(
        "ATOM      1  CA  ALA {}   1       0.000   0.000   0.000  1.00  0.00           C\n",
        chain
    );
    fs::write(&path, line).unwrap();
    path
}

/// A ligand mode shifting every atom 2 A along x, under a 180 degrees rotation around z,
/// must move the ligand 2 A towards -x
macro_rules! anm_rotates_with_ligand {
    ($name:ident, $function:ident, $receptor:expr, $ligand:expr) => {
        #[test]
        fn $name() {
            let mode: Vec<f64> = (0..open($ligand).atom_count())
                .flat_map(|_| vec![2.0, 0.0, 0.0])
                .collect();
            let build = |modes: Vec<f64>, num_anm: usize| {
                $function::new(
                    open($receptor),
                    vec![],
                    vec![],
                    vec![],
                    0,
                    open($ligand),
                    vec![],
                    vec![],
                    modes,
                    num_anm,
                    num_anm > 0,
                )
            };
            let rotation = Quaternion::new(0.0, 0.0, 0.0, 1.0);

            let flexible = build(mode, 1).energy(&[0.0, 0.0, 0.0], &rotation, &[], &[1.0]);
            let expected = build(vec![], 0).energy(&[-2.0, 0.0, 0.0], &rotation, &[], &[]);

            assert!((flexible - expected).abs() < 1e-8);
        }
    };
}

/// Restraints on both partners triple the score only if their atoms are in contact
macro_rules! contact_cutoff {
    ($name:ident, $function:ident) => {
        #[test]
        fn $name() {
            let receptor = single_atom(concat!(stringify!($name), "_receptor"), 'A');
            let ligand = single_atom(concat!(stringify!($name), "_ligand"), 'B');
            let build = |rec_restraints: Vec<String>, lig_restraints: Vec<String>| {
                $function::new(
                    open(&receptor),
                    rec_restraints,
                    vec![],
                    vec![],
                    0,
                    open(&ligand),
                    lig_restraints,
                    vec![],
                    vec![],
                    0,
                    false,
                )
            };
            let free = build(vec![], vec![]);
            let restrained = build(vec!["A.ALA.1".to_string()], vec!["B.ALA.1".to_string()]);

            for (distance, in_contact) in [
                (2.0, true),
                (3.0, true),
                (3.89, true),
                (3.91, false),
                (5.0, false),
            ] {
                let translation = [distance, 0.0, 0.0];
                let rotation = Quaternion::default();
                let score = free.energy(&translation, &rotation, &[], &[]);
                let biased = restrained.energy(&translation, &rotation, &[], &[]);
                let factor = if in_contact { 3.0 } else { 1.0 };

                assert!(score != 0.0, "score is zero at {} A", distance);
                assert!(
                    (biased - factor * score).abs() < 1e-8,
                    "wrong contact at {} A: score {} restrained {}",
                    distance,
                    score,
                    biased
                );
            }
        }
    };
}

anm_rotates_with_ligand!(
    dfire_anm_rotates_with_ligand,
    DFIRE,
    "tests/2oob/2oob_receptor.pdb",
    "tests/2oob/2oob_ligand.pdb"
);
anm_rotates_with_ligand!(
    dna_anm_rotates_with_ligand,
    DNA,
    "tests/1azp/1azp_receptor.pdb",
    "tests/1azp/1azp_ligand.pdb"
);
anm_rotates_with_ligand!(
    pydock_anm_rotates_with_ligand,
    PYDOCK,
    "tests/1azp/1azp_receptor.pdb",
    "tests/1azp/1azp_ligand.pdb"
);

contact_cutoff!(dfire_contact_cutoff, DFIRE);
contact_cutoff!(dna_contact_cutoff, DNA);
contact_cutoff!(pydock_contact_cutoff, PYDOCK);
