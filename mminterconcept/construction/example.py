'''
Author: Andrew-AbiMansour
DOC: Oct 21, 2019

Testing script
'''

import os
from gmx_api import Gmx

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


def solvate_ionize_protein(mdp, salinity, **args):
	""" Returns mdtraj.Trajectory """
	Eng = initialize(mdp, **args)
	Eng = solvate(Eng, mdp)
	Eng = ionize(Eng, salinity, mdp)

	System, top = Eng.getSystem(), Eng.getTop()
	Eng.clean()

	return System, top


def solvate_protein(mdp, **args):
	Eng = initialize(mdp, **args)
	Eng = solvate(Eng, mdp)

	System, top = Eng.getSystem(), Eng.getTop()
	Eng.clean()

	return System, top

def vacuum_protein(mdp, **args):
	Eng = initialize(mdp, **args)

	System, top = Eng.getSystem(), Eng.getTop()
	Eng.clean()

	return System, top

if __name__ == '__main__':
	try:
		pdbID = '2PVP' #'1L2Y'
		mdp = os.path.abspath('../data/1LFH/gmx/mdp/em.mdp')
		pdbFile = os.path.abspath('../data/1LFH/gmx/struct/dialanine.pdb')
		clean = False

		#########################################
		####### Protein in a box of vacuum ######
		#########################################
		System_vacuum, vacuum_top = vacuum_protein(pdbID=pdbID, mdp=mdp, fdir='vacuum', clean=clean)

		########################################
		####### Protein in a box of water ######
		System_solvated, solvated_top = solvate_protein(pdbFile=pdbFile, mdp=mdp, fdir='solvated', clean=clean)

		# P.S. System_ionized without solvent is unphysical / not supported

		########################################
		### Protein in a box of water + ions ###
		########################################
		System_solvated_ionized, ionized_top = solvate_ionize_protein(pdbID=pdbID, mdp=mdp,
			salinity=0.1, fdir='solvated_ionized', clean=clean)

	except Exception:
		raise
