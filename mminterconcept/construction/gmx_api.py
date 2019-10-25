'''
Author: Andrew-AbiMansour
DOC: Oct 23, 2019

GMX (Gromacs) API  for setting up an MD simulation
'''

from api import Workunit, Engine
from collections import OrderedDict
from tools import RandString
import os

class Gmx(Engine):

	def __init__(self, pdbID : str, ff_solute='amber99', ff_solvent='tip3p', topfname='topol.top',
			tprfname='topol.tpr', posrename='posre.itp', ofname='system.gro', ext='gro', 
			exec='gmx', **args):

		args['fdir'] = 'GMX_' + RandString.name()
		super().__init__(pdbID=pdbID, ff_solute=ff_solute, ff_solvent=ff_solvent, topfname=topfname, 
				ofname=ofname, ext=ext, exec=exec, **args)

		self._args['posrename'] = posrename
		self._args['tprfname'] = tprfname
		self._keep = True

	def _updateArgs(self, **rargs):
		for key in rargs:
			if key == 'struct':
				self._args['ifname'] = rargs[key]
			elif key == 'top':
				if isinstance(rargs['top'], list):
					self._args['topfname'] = rargs[key][0]
					self._args['posrename'] = rargs[key][1]
				else:
					self._args['topfname'] = rargs[key]

	def genTop(self):
		with Workunit(keep=self._keep) as Pdb2gmx:
			gmx_args = OrderedDict()

			gmx_args['_exec'] = self._args['_exec']
			gmx_args['_tool'] = 'pdb2gmx'

			gmx_args['-ff'] = self._args['ff_solute']
			gmx_args['-water'] = self._args['ff_solvent']
			gmx_args['-f'] = self._args['ifname']
			gmx_args['-o'] = self._args['ofname']
			gmx_args['-p'] = self._args['topfname']
			gmx_args['-i'] = self._args['posrename']

			args = self._dict_to_str(**gmx_args)

			Pdb2gmx.run(cmd=args)
			rargs = Pdb2gmx.store(sdir=self._sdir, struct=gmx_args['-o'], top=[gmx_args['-p'], gmx_args['-i']])

			self._updateArgs(**rargs)

	def buildBox(self, center=False):
		with Workunit(keep=self._keep) as Editconf:

			gmx_args = OrderedDict()
			gmx_args['_exec'] = self._args['_exec']
			gmx_args['_tool'] = 'editconf'
			gmx_args['-f'] = self._args['ifname']
			gmx_args['-o'] = self._args['ofname']
			gmx_args['-box'] = ('{} ' * len(self._Box.bound)).format(*self._Box.bound)

			if center:
				gmx_args['-c'] = ' '

			args = self._dict_to_str(**gmx_args)

			Editconf.run(cmd=args)
			rargs = Editconf.store(sdir=self._sdir, struct=gmx_args['-o'])

			self._updateArgs(**rargs)

	def solvate(self):
		with Workunit(keep=self._keep) as Solvate:

			gmx_args = OrderedDict()
			gmx_args['_exec'] = self._args['_exec']
			gmx_args['_tool'] = 'solvate'
			gmx_args['-cp'] = self._args['ifname']
			gmx_args['-cs'] = 'spc216'
			gmx_args['-p'] = self._args['topfname']
			gmx_args['-o'] = self._args['ofname']

			args = self._dict_to_str(**gmx_args)

			Solvate.run(cmd=args, input='SOL')
			rargs = Solvate.store(sdir=self._sdir, struct=gmx_args['-o'], top=gmx_args['-p'])

			self._updateArgs(**rargs)

	def genConfig(self, mdpfile, maxwarn=0):

		with Workunit(keep=self._keep) as Grompp:
			gmx_args = OrderedDict()
			gmx_args['_exec'] = self._args['_exec']
			gmx_args['_tool'] = 'grompp'
			gmx_args['-f'] = mdpfile
			gmx_args['-c'] = self._args['ifname']
			gmx_args['-p'] = self._args['topfname']
			gmx_args['-pp'] = 'solvated.top'
			gmx_args['-o'] = self._args['tprfname']
			gmx_args['-maxwarn'] = str(maxwarn)

			print(os.getcwd())
			args = self._dict_to_str(**gmx_args)
			Grompp.run(cmd=args)
			rargs = Grompp.store(sdir=self._sdir, tpr=gmx_args['-o'], top=gmx_args['-pp'])

			self._updateArgs(**rargs)

if __name__ == '__main__':
	Test = Gmx(pdbID='1LFH', box=(9.5, 9, 7))
	Test.genTop()
	Test.buildBox(center=True)
	Test.solvate()
	Test.genConfig(mdpfile='../../../data/1LFH/gmx/mdp/em.mdp', maxwarn=1)
