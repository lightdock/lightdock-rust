use lib3dmol::structures::Structure;
use std::collections::HashMap;
use super::constants::INTERFACE_CUTOFF;


macro_rules! hashmap {
    ($( $key: expr => $val: expr ),*) => {{
         let mut map = ::std::collections::HashMap::new();
         $( map.insert($key, $val); )*
         map
    }}
}


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

pub fn r3_to_numerical(residue_name: &str) -> usize {
    match residue_name {
        "ALA" => {0}
        "CYS" => {1}
        "ASP" => {2}
        "GLU" => {3}
        "PHE" => {4}
        "GLY" => {5}
        "HIS" => {6}
        "ILE" => {7}
        "LYS" => {8}
        "LEU" => {9}
        "MET" => {10}
        "ASN" => {11}
        "PRO" => {12}
        "GLN" => {13}
        "ARG" => {14}
        "SER" => {15}
        "THR" => {16}
        "VAL" => {17}
        "TRP" => {18}
        "TYR" => {19}
        "MMB" => {20}
        _ => { panic!("Residue name not supported in DFIRE scoring function") }
    }
}

// DFIRE only uses 20 distance bins
const DIST_TO_BINS: &[usize] = &[1, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 14, 15, 15, 
                                 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22, 22, 23, 23, 
                                 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31];

lazy_static! {
    static ref ATOMNUMBER: HashMap<&'static str, usize> = hashmap![
        "ALAN" => 0, "ALACA" => 1, "ALAC" => 2, "ALAO" => 3, "ALACB" => 4, 
        "CYSN" => 0, "CYSCA" => 1, "CYSC" => 2, "CYSO" => 3, "CYSCB" => 4, "CYSSG" => 5, 
        "ASPN" => 0, "ASPCA" => 1, "ASPC" => 2, "ASPO" => 3, "ASPCB" => 4, "ASPCG" => 5, "ASPOD1" => 6, "ASPOD2" => 7, 
        "GLUN" => 0, "GLUCA" => 1, "GLUC" => 2, "GLUO" => 3, "GLUCB" => 4, "GLUCG" => 5, "GLUCD" => 6, "GLUOE1" => 7, "GLUOE2" => 8, 
        "PHEN" => 0, "PHECA" => 1, "PHEC" => 2, "PHEO" => 3, "PHECB" => 4, "PHECG" => 5, "PHECD1" => 6, "PHECD2" => 7, "PHECE1" => 8, "PHECE2" => 9, "PHECZ" => 10, 
        "GLYN" => 0, "GLYCA" => 1, "GLYC" => 2, "GLYO" => 3, 
        "HISN" => 0, "HISCA" => 1, "HISC" => 2, "HISO" => 3, "HISCB" => 4, "HISCG" => 5, "HISND1" => 6, "HISCD2" => 7, "HISCE1" => 8, "HISNE2" => 9, 
        "ILEN" => 0, "ILECA" => 1, "ILEC" => 2, "ILEO" => 3, "ILECB" => 4, "ILECG1" => 5, "ILECG2" => 6, "ILECD1" => 7, 
        "LYSN" => 0, "LYSCA" => 1, "LYSC" => 2, "LYSO" => 3, "LYSCB" => 4, "LYSCG" => 5, "LYSCD" => 6, "LYSCE" => 7, "LYSNZ" => 8, 
        "LEUN" => 0, "LEUCA" => 1, "LEUC" => 2, "LEUO" => 3, "LEUCB" => 4, "LEUCG" => 5, "LEUCD1" => 6, "LEUCD2" => 7, 
        "METN" => 0, "METCA" => 1, "METC" => 2, "METO" => 3, "METCB" => 4, "METCG" => 5, "METSD" => 6, "METCE" => 7, 
        "ASNN" => 0, "ASNCA" => 1, "ASNC" => 2, "ASNO" => 3, "ASNCB" => 4, "ASNCG" => 5, "ASNOD1" => 6, "ASNND2" => 7, 
        "PRON" => 0, "PROCA" => 1, "PROC" => 2, "PROO" => 3, "PROCB" => 4, "PROCG" => 5, "PROCD" => 6, 
        "GLNN" => 0, "GLNCA" => 1, "GLNC" => 2, "GLNO" => 3, "GLNCB" => 4, "GLNCG" => 5, "GLNCD" => 6, "GLNOE1" => 7, "GLNNE2" => 8, 
        "ARGN" => 0, "ARGCA" => 1, "ARGC" => 2, "ARGO" => 3, "ARGCB" => 4, "ARGCG" => 5, "ARGCD" => 6, "ARGNE" => 7, "ARGCZ" => 8, "ARGNH1" => 9, "ARGNH2" => 10, 
        "SERN" => 0, "SERCA" => 1, "SERC" => 2, "SERO" => 3, "SERCB" => 4, "SEROG" => 5, 
        "THRN" => 0, "THRCA" => 1, "THRC" => 2, "THRO" => 3, "THRCB" => 4, "THROG1" => 5, "THRCG2" => 6, 
        "VALN" => 0, "VALCA" => 1, "VALC" => 2, "VALO" => 3, "VALCB" => 4, "VALCG1" => 5, "VALCG2" => 6, 
        "TRPN" => 0, "TRPCA" => 1, "TRPC" => 2, "TRPO" => 3, "TRPCB" => 4, "TRPCG" => 5, "TRPCD1" => 6, "TRPCD2" => 7, "TRPCE2" => 8, "TRPNE1" => 9, "TRPCE3" => 10, "TRPCZ3" => 11, "TRPCH2" => 12, "TRPCZ2" => 13, 
        "TYRN" => 0, "TYRCA" => 1, "TYRC" => 2, "TYRO" => 3, "TYRCB" => 4, "TYRCG" => 5, "TYRCD1" => 6, "TYRCD2" => 7, "TYRCE1" => 8, "TYRCE2" => 9, "TYRCZ" => 10, "TYROH" => 11, 
        "MMBBJ" => 0];

    // Atom type and residue translation matrix
    static ref ATOMRES: Vec<Vec<usize>> = vec![vec![74, 75, 76, 77, 78, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![0, 1, 2, 3, 4, 5, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![122, 123, 124, 125, 126, 127, 128, 129, 0, 0, 0, 0, 0, 0],
                                               vec![113, 114, 115, 116, 117, 118, 119, 120, 121, 0, 0, 0, 0, 0],
                                               vec![14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 0, 0, 0],
                                               vec![79, 80, 81, 82, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 0, 0, 0, 0],
                                               vec![25, 26, 27, 28, 29, 30, 31, 32, 0, 0, 0, 0, 0, 0],
                                               vec![151, 152, 153, 154, 155, 156, 157, 158, 159, 0, 0, 0, 0, 0],
                                               vec![33, 34, 35, 36, 37, 38, 39, 40, 0, 0, 0, 0, 0, 0],
                                               vec![6, 7, 8, 9, 10, 11, 12, 13, 0, 0, 0, 0, 0, 0],
                                               vec![105, 106, 107, 108, 109, 110, 111, 112, 0, 0, 0, 0, 0, 0],
                                               vec![160, 161, 162, 163, 164, 165, 166, 0, 0, 0, 0, 0, 0, 0],
                                               vec![96, 97, 98, 99, 100, 101, 102, 103, 104, 0, 0, 0, 0, 0],
                                               vec![140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 0, 0, 0],
                                               vec![90, 91, 92, 93, 94, 95, 0, 0, 0, 0, 0, 0, 0, 0],
                                               vec![83, 84, 85, 86, 87, 88, 89, 0, 0, 0, 0, 0, 0, 0],
                                               vec![41, 42, 43, 44, 45, 46, 47, 0, 0, 0, 0, 0, 0, 0],
                                               vec![48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61],
                                               vec![62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 0, 0],
                                               vec![167, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]];
}


#[derive(Debug)]
pub struct DockingModel {
    pub atoms: Vec<usize>,
    pub coordinates: Vec<[f64; 3]>,
    pub membrane: Vec<usize>,
    pub active_restraints: HashMap<String, Vec<usize>>, 
    pub passive_restraints: HashMap<String, Vec<usize>>,
    pub nmodes: Vec<f64>,
    pub num_anm: usize,
}


impl DockingModel {
    pub fn new(structure: &Structure, active_restraints: &[String], passive_restraints: &[String],
        nmodes: &Vec<f64>, num_anm: usize) -> Self {
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
                let res_id = format!("{}.{}.{}", chain.name, residue.name.trim(), residue.res_num);

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
                    
                    let rnuma = r3_to_numerical(residue.name.trim());
                    let anuma = match ATOMNUMBER.get(&rec_atom_type[..]) {
                        Some(&a) => a as usize,
                        _ => panic!("Not supported atom type {:?}", rec_atom_type),
                    };
                    let atoma = ATOMRES[rnuma][anuma];
                    model.atoms.push(atoma);
                    model.coordinates.push([atom.coord[0] as f64, atom.coord[1] as f64, atom.coord[2] as f64]);
                    atom_index += 1;
                }
            }
        }
        model
    }
}


#[derive(Debug)]
pub struct DFIRE {
    pub energy: Vec<f64>,
}

impl Default for DFIRE {
    fn default() -> Self {
        DFIRE::new()
    }
}


impl DFIRE {
    pub fn new() -> DFIRE {
        let mut d = DFIRE {
            energy: Vec::with_capacity(168 * 168 * 20)
        };
        d.load_potentials();
        d
    }

    fn load_potentials(&mut self) {
        let raw_parameters = include_str!("DCparams").to_string();
        let split = raw_parameters.lines();
        let params: Vec<&str> = split.collect();

        for count in 0..168*168*20 {
            self.energy.push(params[count].trim().parse::<f64>().unwrap());
        }
    }

    pub fn get_potential(&mut self, x: usize, y: usize, z: usize) -> f64{
        self.energy[x + 168 * (y + 20 * z)]
    }

    pub fn get_docking_model(structure: &Structure, 
        active_restraints: &[String], passive_restraints: &[String], 
        nmodes: &Vec<f64>, num_anm: usize) -> DockingModel {
        DockingModel::new(&structure, &active_restraints, &passive_restraints, &nmodes, num_anm)
    }

    pub fn energy(&self, receptor: &DockingModel, ligand: &DockingModel,
        receptor_coordinates: &[[f64; 3]], ligand_coordinates: &[[f64; 3]],
        interface_receptor: &mut Vec<usize>, interface_ligand: &mut Vec<usize>) -> f64 {
        let mut score: f64 = 0.0;
        for (i, ra) in receptor_coordinates.iter().enumerate() {
            let x1 = ra[0];
            let y1 = ra[1];
            let z1 = ra[2];
            let atoma = receptor.atoms[i];
            for (j, la) in ligand_coordinates.iter().enumerate() {
                let dist = (x1-la[0])*(x1-la[0]) + (y1-la[1])*(y1-la[1]) + (z1-la[2])*(z1-la[2]);
                if dist <= 225. {
                    let atomb = ligand.atoms[j];
                    let d = dist.sqrt()*2.0 - 1.0;
                    let dfire_bin = DIST_TO_BINS[d as usize] - 1;
                    score += self.energy[atoma*168*20 + atomb*20 + dfire_bin];
                    if d <= INTERFACE_CUTOFF {
                        interface_receptor[i] = 1;
                        interface_ligand[j] = 1;
                    }
                }
            }
        }
        (score * 0.0157 - 4.7) * -1.0
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_read_potentials() {
        let scoring = DFIRE::new();
        assert_eq!(scoring.energy[0], 10.0);
        assert_eq!(scoring.energy[2], -0.624030868);
        assert_eq!(scoring.energy[4998], -0.0458685914);
        assert_eq!(scoring.energy[168*168*20-1], 0.0);
    }
}
