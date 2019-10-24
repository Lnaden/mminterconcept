'''
Author: Andrew-AbiMansour
DOC: Oct 23, 2019

GMX (Gromacs) API  for setting up an MD simulation
'''

from api import Workunit, Engine
from collections import OrderedDict
import os

class Gmx(Engine):

	def __init__(self, pdbID : str, ff_solute='amber99', ff_solvent='tip3p', topfname='topol.top',
			posrename='posre.itp', ofname='system.gro', ext='gro', **args):

		super().__init__(pdbID=pdbID, ff_solute=ff_solute, ff_solvent=ff_solvent, topfname=topfname, 
				ofname=ofname, ext=ext, **args)

		self._gmx_args = {
			'_exec': 'gmx',
			'ifname': os.path.abspath(f'{pdbID}.{ext}'),
			'ff_solute': ff_solute,
			'ff_solvent': ff_solvent,
			'ofname': ofname,
			'topfname': topfname,
			'posrename': posrename
		}

	def pdb2gmx(self):
		with Workunit(keep=False) as Pdb2gmx:
			gmx_args = OrderedDict()

			# HACKISH!!
			gmx_args['_exec'] = self._gmx_args['_exec']
			gmx_args['_tool'] = 'pdb2gmx'

			gmx_args['-ff'] = self._gmx_args['ff_solute']
			gmx_args['-water'] = self._gmx_args['ff_solvent']
			gmx_args['-f'] = self._gmx_args['ifname']
			gmx_args['-o'] = self._gmx_args['ofname']
			gmx_args['-p'] = self._gmx_args['topfname']
			gmx_args['-i'] = self._gmx_args['posrename']

			args = self._dict_to_str(**gmx_args)

			Pdb2gmx.run(cmd=args)
			Pdb2gmx.store(sdir=self._sdir, struct=gmx_args['-o'], top=[gmx_args['-p'], gmx_args['-i']])

			self._gmx_args['ifname'] =  os.path.abspath(os.path.join(self._sdir, 'struct', gmx_args['-o']))
			self._gmx_args['topfname'] =  os.path.abspath(os.path.join(self._sdir, 'top', gmx_args['-p']))
			self._gmx_args['posrename'] =  os.path.abspath(os.path.join(self._sdir, 'top', gmx_args['-i']))

	def editconf(self):
		with Workunit(keep=False) as Editconf:

			gmx_args = OrderedDict()
			gmx_args['_exec'] = self._gmx_args['_exec']
			gmx_args['_tool'] = 'editconf'
			gmx_args['-f'] = self._gmx_args['ifname']
			gmx_args['-o'] = self._gmx_args['ofname']
			gmx_args['-c'] = ' '
			gmx_args['-box'] = ('{} ' * len(self._Box.bound)).format(*self._Box.bound)

			args = self._dict_to_str(**gmx_args)
			print(args)

			Editconf.run(cmd=args)
			Editconf.store(sdir=self._sdir, struct=gmx_args['-o'])

			self._gmx_args['ifname'] = os.path.abspath(os.path.join(self._sdir, 'struct', gmx_args['-o']))

	def solvate(self):
		with Workunit(keep=False) as Solvate:

			gmx_args = OrderedDict()
			gmx_args['_exec'] = self._gmx_args['_exec']
			gmx_args['_tool'] = 'solvate'
			gmx_args['-cp'] = self._gmx_args['ifname']
			gmx_args['-cs'] = 'spc216'
			gmx_args['-p'] = self._gmx_args['topfname']
			gmx_args['-o'] = self._gmx_args['ofname']

			args = self._dict_to_str(**gmx_args)
			print(args)

			Solvate.run(cmd=args, input='SOL')
			Solvate.store(sdir=self._sdir, struct=gmx_args['-o'], top=gmx_args['-p'])

			self._gmx_args['ifname'] =  os.path.abspath(os.path.join(self._sdir, 'struct', gmx_args['-o']))
			self._gmx_args['topfname'] =  os.path.abspath(os.path.join(self._sdir, 'top', gmx_args['-p']))
			self._gmx_args['posrename'] =  os.path.abspath(os.path.join(self._sdir, 'top'))

if __name__ == '__main__':
	Test = Gmx(pdbID='1LFH', box=(9.5, 9, 7))
	Test.pdb2gmx()
	Test.editconf()
	Test.solvate()
