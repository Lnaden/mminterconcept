'''
Author: Andrew-AbiMansour
DOC: Oct 21, 2019

Testing script (eventually will be deleted)
'''

import sys
from gmx_api import Gmx

def run(pdbID, mdp='../../../data/1LFH/gmx/mdp/em.mdp'):

	Test = Gmx(pdbID=pdbID, box=(9.5, 9, 7))

	# Run pdb2gmx to generate topology
	Test.genTop()

	# Construct sim box in which the protein is centered 
	Test.buildBox(center=True)

	# Solvate box with water using the spc216 model
	Test.solvate(model='spc216')

	# Generate binary (tpr) config file from em.mdp
	Test.genConfig(mdpfile=mdp, maxwarn=1)

	# Add NaCL (default) ions of salinity 0.1M
	Test.addIons(conc=0.1, replace='SOL')

	# Update tpr config file
	Test.genConfig(mdpfile=mdp)

	# Return mdtraj.Trajectory
	return Test.load()

if __name__ == '__main__':
	try:
		pdbID = sys.argv[1]
		System = run(pdbID)
	except Exception:
		raise
