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

class ImportPDB:

        def __init__(self, pdbID, length=6):
                self.fname = ''.join(random.SystemRandom().choice(string.ascii_uppercase + string.digits) for _ in range(length)) + '.pdb'
                self.pdbID = pdbID

        def __enter__(self):
                with open(self.fname,'w') as fp:
                        fp.write(pypdb.get_pdb_file(self.pdbID))

                return mdtraj.load(self.fname)

        def __exit__(self, exc_type, exc_val, exc_tb):
                os.remove(self.fname)
