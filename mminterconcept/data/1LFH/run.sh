#!/usr/bin/env bash
# requirements: Gromacs 5.0 and above

export DIR=.
export EXT='pdb'

rm $DIR/struct/*
rm $DIR/top/*

# Convert pdb to gro using amber99 FF for 1LFH and TIP3P FF for water
gmx pdb2gmx -f ../1LFH.pdb -ff amber99 -water tip3p -o $DIR/struct/protein.${EXT} -p $DIR/top/topol.top -i $DIR/top/posre.itp

# Center protein in simulation box of size (9.5,9,7)nm
gmx editconf -f $DIR/struct/protein.${EXT} -c -o $DIR/struct/protein.${EXT} -box 9.5 9 7

# Solvate system based on the gromacs-SPC model
gmx solvate -cp $DIR/struct/protein.${EXT} -cs spc216 -p $DIR/top/topol.top -o $DIR/struct/protein_solvated.${EXT}

# Generate binary (tpr) output file from input data and EM mdp file
gmx grompp -f $DIR/mdp/em.mdp -c $DIR/struct/protein_solvated.${EXT} -p $DIR/top/topol.top -pp $DIR/top/system_solvated.top -po $DIR/mdp/em.mdp  -o $DIR/topol.tpr -maxwarn 1

# Add ions
gmx genion -s $DIR/topol.tpr -p $DIR/top/system_solvated.top -conc 0.1 -o $DIR/struct/system_ionized.${EXT}  -neutral

# Update tpr file
gmx grompp -f $DIR/mdp/em.mdp -c $DIR/struct/system_ionized.${EXT} -p $DIR/top/system_solvated.top -pp $DIR/top/system_ionized.top -po $DIR/mdp/em.mdp  -o $DIR/topol.tpr

