'''
Author: Andrew-AbiMansour
DOC: Oct 22, 2019

Generic API for MD solvers (e.g. Gromacs)
'''

import subprocess
from tools import RandString, ImportPDB
from collections import OrderedDict
import os, sys
import errno
import shutil

class Engine:
	def __init__(self, pdbID : str, ff_solute: str, ff_solvent: str, topfname: str, ofname: str, ext: str, **args):

		if 'fdir' in args:
			self._fdir = args['fdir']
		else:
			self._fdir = 'GMX_' + RandString.name()

		if os.path.isdir(self._fdir):
			raise Exception(f'__enter__ : [ERROR]: randrom dir {self._fdir} already exists')
		else:
			os.mkdir(self._fdir)
			os.chdir(self._fdir)

		if 'sdir' in args:
			self._sdir = sdir
		else:
			self._sdir = '..'

		self._pdbID = pdbID

		# dict not always ordered (depends on python ver)
		self._gmx_args = OrderedDict()

		with ImportPDB(self._pdbID) as SS:
			SS.save(f'{pdbID}.{ext}')

		if not os.path.exists('struct'):
			os.mkdir('struct')

		if not os.path.exists('top'):
			os.mkdir('top')

		if not os.path.exists('top'):
			os.mkdir('mdp')

	def _dict_to_flist(self, **args) -> list:
		""" Map dict to flattened list """
		args = [[key, val] for key, val in args.items()]
		args = [single for pair in args for single in pair]
		return [arg for arg in args if not arg.startswith('_')]

class Workunit:
	def __init__(self, keep=False, fdir=None):
		self._keep = keep
		if fdir:
			self._fdir = fdir
		else:
			self._fdir = 'Unit_' + RandString.name()

	def __enter__(self):
		if os.path.isdir(self._fdir):
			raise Exception(f'__enter__ : [ERROR]: randrom dir {self._fdir} already exists')
		else:
			os.mkdir(self._fdir)
			os.chdir(self._fdir)

		return self

	def run(self, cmd, wait=False):
		try:
			if (wait):
				proc = subprocess.Popen(cmd, stdout = subprocess.PIPE)
				proc.wait()
			else:
				proc = subprocess.Popen(cmd, stdin = None, stdout = None, stderr = None, close_fds = True)

			(result, error) = proc.communicate()

		except subprocess.CalledProcessError as err:
			sys.stderr.write(f'run : [ERROR]: output = {err.output}, error code = {err.returncode}')

		return result

	def __exit__(self, exc_type, exc_val, exc_tb):
		os.chdir('..')
		self.clean()

	def store(self, sdir, **files):
		try:
			for fdir, sfile in files.items():
				if not os.path.exists(fdir):
					os.mkdir(fdir)

				if isinstance(sfile, str):
					if(os.path.isfile(sfile)):
						if os.path.exists(os.path.join(sdir, fdir)):
							shutil.move(sfile, os.path.join(sdir, fdir, sfile))
						else:
							raise 

				elif isinstance(sfile, list):
					for ssfile in sfile:
						if(os.path.isfile(ssfile)):
							if os.path.exists(os.path.join(sdir, fdir)):
								shutil.move(ssfile, os.path.join(sdir, fdir, ssfile))
							else:
								raise
						else:
							raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), ssfile)
				else:
					raise TypeError(f'store : [ERROR]: sfile = {type(sfile)}, must be str or list')

		except Exception:
			raise

	def clean(self):
		if not self._keep:
			if os.path.isdir(self._fdir):
				shutil.rmtree(self._fdir, ignore_errors=True)
