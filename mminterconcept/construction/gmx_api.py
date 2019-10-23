'''
Author: Andrew-AbiMansour
DOC: Oct 22, 2019

Generic API for MD solvers (e.g. Gromacs)
'''

from api import Process
from tools import ImportPDB

with Process(keep=True) as Pdb2gmx:

	speciesfname = '1LFH.gro'
	species_ff = 'amber99'
	water_ff = 'tip3p'

	with ImportPDB(pdbID='1LFH') as SS:
		SS.save(speciesfname)

	Pdb2gmx.run(cmd=['gmx', 'pdb2gmx', '-f', f'{speciesfname}', '-ff', species_ff, '-water', water_ff, '-o', 'protein.gro', '-p', 'topol.top', '-i', 'posre.itp'])
	Pdb2gmx.store(struct='protein.gro', top=['topol.top','posre.itp'])
