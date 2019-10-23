'''
Author: Andrew-AbiMansour
DOC: Oct 21, 2019

Testing script (eventually will be deleted)
'''

import sys
from tools import ImportPDB
import core

def run(pdbID):
	with ImportPDB(pdbID) as Spec:
        	Species = Spec

	Box = core.Box(bound = ( (0,1),(0,1),(0,1) ))
	Sys = core.System(Box=Box, Species=Species)

if __name__ == '__main__':
	try:
		pdbID = sys.argv[1]
		run(pdbID)
	except Exception:
		raise
