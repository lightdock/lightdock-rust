#!/bin/bash

# Edit only the number of cores to use for this simulation
NUM_CORES=4

# LightDock setup
lightdock3_setup.py 1czy_protein.pdb 1czy_peptide.pdb 400 200 --noxt --noh -anm -rst restraints.list

# Convert ANM data
lgd_flatten.py lightdock_rec.nm.npy rec_nm.npy
lgd_flatten.py lightdock_lig.nm.npy lig_nm.npy

# Calculate number of swarms
s=`ls -d swarm_* | wc -l`
swarms=$((s-1))

# Copy binary
cp ../../target/release/lightdock-rust .

# Create a task.list file for ant_thony
for i in `seq 0 $swarms`;do echo "cd swarm_${i}; cp ../lightdock_1czy_protein.pdb .; cp ../lightdock_1czy_peptide.pdb .;cp ../rec_nm.npy .;cp ../lig_nm.npy .;cp -R ../../../data .;../lightdock-rust ../setup.json ../init/initial_positions_${i}.dat 10 dfire; rm -rf lightdock_*.pdb *.npy data;" >> task.list; done

# Let ant_thony run
ant_thony.py --cores ${NUM_CORES} task.list

# Clean task.list
rm -rf task.list
