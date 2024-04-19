# LightDock-Rust

A Rust implementation of the [LightDock](https://lightdock.org) macromolecular software with the DFIRE and DNA scoring functions.

## Installation
1. Clone this repository:

 ```
 git clone https://github.com/lightdock/lightdock-rust.git
 ```

2. Compile it with Rust (you may install Rust using [rustup](https://rustup.rs/)):

 ```
 cd lightdock-rust
 cargo build --release
 ```
 
## Examples

Several examples can be found in the `example` folder.

| :information_source: Data Path          |
|:---------------------------|
| You may set an environment variable `LIGHTDOCK_DATA` to point to the data folder included in this repository to avoid copying it: `export LIGHTDOCK_DATA=/path/to/lightdock-rust/data`  |

Recorded times on MacBook Pro M3 Pro.

### 1k4c (Membrane docking)

```bash
cd example/1k4c
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire
```

Output:

```
Reading starting positions from "initial_positions_0.dat"
Swarm ID 0
Writing to swarm dir "swarm_0"
Reading receptor input structure: lightdock_receptor_membrane.pdb
Reading ligand input structure: lightdock_ligand.pdb
Loading DFIRE scoring function
Creating GSO with 200 glowworms
Starting optimization (100 steps)

real    1m52,132s
user    1m51,808s
sys     0m0,150s
```

### 1ppe (protein docking)

```bash
cd example/1ppe
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire
```

Output:

```
Reading starting positions from "initial_positions_0.dat"
Swarm ID 0
Writing to swarm dir "swarm_0"
Reading receptor input structure: lightdock_1ppe_e.pdb
Reading ligand input structure: lightdock_1ppe_i.pdb
Loading DFIRE scoring function
Creating GSO with 200 glowworms
Starting optimization (100 steps)

real    0m4,252s
user    0m4,142s
sys     0m0,092s
```

### 2uuy (protein docking)

```bash
cd example/2uuy
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire
```

Output:

```
Reading starting positions from "initial_positions_0.dat"
Swarm ID 0
Output directory does not exist for swarm 0, creating it...
Writing to swarm dir "swarm_0"
Reading receptor input structure: lightdock_2UUY_rec.pdb
Reading ligand input structure: lightdock_2UUY_lig.pdb
Loading DFIRE scoring function
Creating GSO with 200 glowworms
Starting optimization (100 steps)

real    0m8,108s
user    0m7,834s
sys	    0m0,261s
```

### 1czy (protein-peptide docking)

```bash
cd example/1czy
time ../../target/release/lightdock-rust setup.json init/initial_positions_0.dat 100 dfire
```

Output:

```
Reading starting positions from "init/initial_positions_0.dat"
Swarm ID 0
Writing to swarm dir "swarm_0"
Reading receptor input structure: lightdock_1czy_protein.pdb
Reading ligand input structure: lightdock_1czy_peptide.pdb
Loading DFIRE scoring function
Creating GSO with 200 glowworms
Starting optimization (100 steps)

real    0m1,580s
user    0m1,312s
sys	    0m0,248s
```

### 1azp (protein-nucleic docking)

```bash
cd example/1azp
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dna
```

Output:

```
Reading starting positions from "initial_positions_0.dat"
Swarm ID 0
Output directory does not exist for swarm 0, creating it...
Writing to swarm dir "swarm_0"
Reading receptor input structure: lightdock_protein.pdb
Reading ligand input structure: lightdock_dna.pdb
Loading DNA scoring function
Creating GSO with 200 glowworms
Starting optimization (100 steps)

real    0m14,228s
user    0m13,932s
sys     0m0,281s
```

