#!/bin/bash

# Edit only the number of cores to use for this simulation
NUM_CORES=4

# LightDock setup
lightdock3_setup.py 1czy_protein.pdb 1czy_peptide.pdb --noxt --noh -anm -rst restraints.list -spr 10

# Convert ANM data
lgd_flatten.py lightdock_rec.nm.npy rec_nm.npy
lgd_flatten.py lightdock_lig.nm.npy lig_nm.npy

# Calculate number of swarms
s=`ls -d swarm_* | wc -l`
swarms=$((s-1))

# Copy lightdock-rust binary
cp ../../target/release/lightdock-rust .

# Create a task.list file for ant_thony
for i in `seq 0 $swarms`;do echo "./lightdock-rust setup.json init/initial_positions_${i}.dat 100 dfire;" >> task.list; done

# Let ant_thony run
time ant_thony.py --cores ${NUM_CORES} task.list

# Clean task.list
rm -rf task.list
