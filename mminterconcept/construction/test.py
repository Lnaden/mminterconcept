'''
Author: Andrew-AbiMansour
DOC: Oct 21, 2019

Testing script (eventually will be deleted)
'''

import sys
from tools import ImportPDB
import system

def run(pdbID):
	with ImportPDB(pdbID) as Spec:
        	species = Spec

	print(species.xyz)

	box = system.Box(bound = ( (0,1),(0,1),(0,1) ))
	print(box.shape, box.bound)

if __name__ == '__main__':
	try:
		pdbID = sys.argv[1]
		run(pdbID)
	except Exception:
		raise
