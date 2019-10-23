'''
Author: Andrew-AbiMansour
DOC: Oct 23, 2019

GMX (Gromacs) API  for setting up an MD simulation
'''

from api import Process
from tools import ImportPDB
from collections import OrderedDict

with Process(keep=True) as Pdb2gmx:

	# dict not always ordered (depends on python ver)
	gmx_args = OrderedDict()

	# HACKISH!!

	gmx_args['_exec'] = 'gmx'
	gmx_args['_tool'] = 'pdb2gmx'
	gmx_args['-f'] = '1LFH.gro'
	gmx_args['-ff'] = 'amber99'
	gmx_args['-water'] = 'tip3p'
	gmx_args['-o'] = 'protein.gro'
	gmx_args['-p'] = 'topol.top'
	gmx_args['-i'] = 'posre.itp'

	with ImportPDB(pdbID='1LFH') as SS:
		SS.save(gmx_args['-f'])

	args = [[key, val] for key, val in gmx_args.items()]
	args = [single for pair in args for single in pair]

	args = [arg for arg in args if not arg.startswith('_')]

	print(args)

	Pdb2gmx.run(cmd=args)
	Pdb2gmx.store(struct=gmx_args['-f'], top=[gmx_args['-p'], gmx_args['-i']])
