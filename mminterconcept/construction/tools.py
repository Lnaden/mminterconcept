'''
Author: Andrew-AbiMansour
DOC: Oct 21, 2019

This module provides common classes/functions used in other modules
 '''

import pypdb
import mdtraj
import random
import string
import os

class RandString:
	@staticmethod
	def name(length=6):
		return ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(length))

class ImportPDB:

	def __init__(self, pdbID):
		self._fname = RandString.name() + '.pdb'
		self._pdbID = pdbID

	def __enter__(self) -> mdtraj.Trajectory:
		with open(self._fname,'w') as fp:
			fp.write(pypdb.get_pdb_file(self._pdbID))

		return mdtraj.load(self._fname)

	def __exit__(self, exc_type, exc_val, exc_tb):
		os.remove(self._fname)
