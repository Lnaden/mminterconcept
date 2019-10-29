'''
Author: Andrew-AbiMansour
DOC: Oct 21, 2019

Testing script
'''

import os
from .gmx_api import Gmx

def ionize(Eng: Gmx, salinity: float, mdp: str, group: str='SOL') -> Gmx:

	# Generate binary (tpr) config file from mdp file
	Eng.genConfig(mdpfile=mdp, maxwarn=1)

	# Add NaCL (default) ions, salinity in mol/L
	Eng.addIons(salinity, replace=group)

	# Update tpr config file
	Eng.genConfig(mdpfile=mdp)

	return Eng

def solvate(Eng: Gmx, mdp: str) -> Gmx:
	# Solvate box with water using the spc216 model
	Eng.solvate(model='spc216')

	Eng.genConfig(mdpfile=mdp)

	return Eng

def initialize(mdp, **args) -> Gmx:

	Eng = Gmx(mdp=mdp, **args)

	# Run pdb2gmx to generate topology
	Eng.genTop()

	# Construct sim box in which the protein is centered
	Eng.buildBox(center=True)

	# Generate binary tpr config file
	Eng.genConfig(mdpfile=mdp)

	return Eng

#####################################################
############ Examples ###############################


def solvate_ionize_protein(pdbID, mdp, salinity, **args):
	""" Returns mdtraj.Trajectory """
	Eng = initialize(mdp, pdbID=pdbID, **args)
	Eng = solvate(Eng, mdp)
	Eng = ionize(Eng, salinity, mdp)

	return Eng.getSystem(), Eng.getTop()


def solvate_protein(pdbID, mdp, **args):
	Eng = initialize(mdp, pdbID=pdbID, **args)
	Eng = solvate(Eng, mdp)

	return Eng.getSystem(), Eng.getTop()


def vacuum_protein(pdbID, mdp, **args):
	Eng = initialize(mdp, pdbID=pdbID, **args)
	return Eng.getSystem(), Eng.getTop()


if __name__ == '__main__':
	try:
		pdbID = '1LFH'
		mdp = os.path.abspath('../data/em.mdp')

		#########################################
		####### Protein in a box of vacuum ######
		#########################################
		System_vacuum, vacuum_top = vacuum_protein(pdbID, mdp, fdir='vacuum')

		# Expand System sim box and then solvate
		System_vacuum

		########################################
		####### Protein in a box of water ######
		System_solvated, solvated_top = solvate_protein(pdbID, mdp, fdir='solvated')

		# P.S. System_ionized without solvent is unphysical / not supported

		########################################
		### Protein in a box of water + ions ###
		########################################
		System_solvated_ionized, ionized_top = solvate_ionize_protein(pdbID, mdp,
			salinity=0.1, fdir='solvated_ionized')

	except Exception:
		raise
