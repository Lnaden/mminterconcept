'''
Author: Andrew-AbiMansour
DOC: Oct 23, 2019

GMX (Gromacs) API  for setting up an MD simulation
'''

from api import Unit, Engine
import os

class Gmx(Engine):

	def __init__(self, pdbID : str, ff_solute='amber99', ff_solvent='tip3p', topfname='topol.top', 
			posrename='posre.itp', ofname='system.gro', ext='gro', **args):

		super().__init__(pdbID=pdbID, ff_solute=ff_solute, ff_solvent=ff_solvent, topfname=topfname, 
				ofname=ofname, ext=ext, **args)

		self._gmx_args['_exec'] = 'gmx'
		self._gmx_args['_tool'] = None
		self._gmx_args['-f'] = os.path.abspath(f'{pdbID}.{ext}')
		self._gmx_args['-ff'] = ff_solute
		self._gmx_args['-water'] = ff_solvent
		self._gmx_args['-o'] = ofname
		self._gmx_args['-p'] = topfname
		self._gmx_args['-i'] = posrename


	def _dict_to_flist(self, **gmx_args) -> list:
		""" Map dict to flattened list """
		args = [[key, val] for key, val in gmx_args.items()]
		args = [single for pair in args for single in pair]
		return [arg for arg in args if not arg.startswith('_')]

	def pdb2gmx(self):
		with Unit(keep=True) as Pdb2gmx:
			gmx_args = self._gmx_args.copy()
			# HACKISH!!
			gmx_args['_tool'] = 'pdb2gmx'

			args = self._dict_to_flist(**gmx_args)

			Pdb2gmx.run(cmd=args)
			Pdb2gmx.store(sdir='..', struct=gmx_args['-o'], top=[gmx_args['-p'], gmx_args['-i']])

if __name__ == '__main__':
	Test = Gmx(pdbID='1LFH')
	Test.pdb2gmx()
