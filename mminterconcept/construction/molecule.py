'''
Author: Andrew-AbiMansour
DOC: Oct 21, 2019

This module defines a molecule for MD simulation/analysis
'''

from qcelemental.models import molecule as qc_molecule
from pydantic import Schema
from typing import Optional, Tuple, List

class Molecule(qc_molecule.Molecule):

	angles: Optional[List[Tuple[int, int, int]]] = Schema(None,
        	description="The 3-body angle defined between 3 connected atoms. Each entry in this "
        		"list is a Tuple of ``(atom_index_A, atom_index_B, atom_index_C)`` where the ``atom_index`` "
       			"matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``.")

	dihedrals: Optional[List[Tuple[int, int, int]]] = Schema(None,
       		description="The 5-body dihedral angle between 2 bonds originating from different atoms. Each entry in this "
        	"list is a Tuple of ``(atom_index_A, atom_index_B, atom_index_C, ...)`` where the ``atom_index`` "
        	"matches the 0-indexed indices of all other per-atom settings like ``symbols`` and ``real``.")



