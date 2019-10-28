'''
Author: Andrew-AbiMansour
DOC: Oct 22, 2019

Generic API for MD solvers (e.g. Gromacs)
'''

import subprocess
from .tools import RandString, ImportPDB
from .core import Box
import os, sys
import errno
import shlex
import shutil
import mdtraj


class PopenWithInput(subprocess.Popen):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def communicate(self, input=None):
		if input:
			return super().communicate(input.encode('utf-8'))
		else:
			return super().communicate()


class Engine:

	def _abspath(self, *args) -> str:
		return os.path.abspath(os.path.join(self._adir, *args))

	def __init__(self, ff_solute: str, ff_solvent: str, topfname: str, ofname: str, ext: str, exec: str, **args):

		if 'wdir' in args:
			os.chdir(args['wdir'])

		if 'fdir' in args:
			self._fdir = args['fdir']
		else:
			self._fdir = RandString.name()

		if 'mod' not in args:
			self._mod = 0o777

		os.makedirs(self._fdir, exist_ok=True)
		os.chmod(self._fdir, self._mod)

		self._adir = os.path.abspath(self._fdir)

		# os.chdir(self._fdir)

		self._args = {
			'_exec': exec,
			'_system': None,
			'ff_solute': ff_solute,
			'ff_solvent': ff_solvent,
			'ofname': ofname,
			'topfname': topfname
		}

		if 'sdir' in args:
			self._sdir = args["sdir"]
		else:
			self._sdir = '..'

		if 'box' not in args:
			self._Box = Box(shape='cubic', bound=(10.0, 10.0, 10.0))  # we could estimate box lengths from protein size
		else:
			self._Box = Box(shape='cubic', bound=args['box'])

		if 'pdbID' in args:
			with ImportPDB(args['pdbID']) as SS:
				# clean system of any water molecules
				drySS = self.getSystemDry(SS)
				drySS.save(self._abspath(f'{args["pdbID"]}.{ext}'))

			self._args['ifname'] = self._abspath(f'{args["pdbID"]}.{ext}')

		else:
			if 'System' in args:
				args['System'].save(self._abspath(f'System.{ext}'))

				self._args['ifname'] = self._abspath(f'{args["System"]}.{ext}')
			else:
				raise ValueError('pdbID or system keywords must be supplied')


		# create new dirs if requested
		if '_static_dirs' in args:
			for path in args['_static_dirs']:
				path = self._abspath(path)
				if not os.path.exists(path):
					os.mkdir(path)
					os.chmod(path, self._mod)

	def _dict_to_str(self, **args):
		""" Map dict to flattened str """
		args = [[key, val] for key, val in args.items()]
		args = [single for pair in args for single in pair]
		return ' '.join([arg for arg in args if not arg.startswith('_')])

	def getSystem(self) -> mdtraj.Trajectory:
		return mdtraj.load(self._abspath(self._args['_system']))

	def getSystemDry(self, System: mdtraj.Trajectory=None) -> mdtraj.Trajectory:
		sel = 'protein' #'not resname HOH and not resname SOL'

		if System:
			return System.atom_slice(System.topology.select(sel))
		else:
			System = mdtraj.load(self._abspath(self._args['_system']))
			return System.atom_slice(System.topology.select(sel))

class Workunit:
	def __init__(self, keep=False, fdir=None, mod=None):
		self._keep = keep
		if fdir:
			self._fdir = fdir
		else:
			self._fdir = 'Unit_' + RandString.name()

		if not mod:
			self._mod = 0o777

	def __enter__(self):
		if os.path.isdir(self._fdir):
			raise Exception(f'__enter__ : [ERROR]: randrom dir {self._fdir} already exists')
		else:
			os.mkdir(self._fdir)
			os.chmod(self._fdir, self._mod)
			os.chdir(self._fdir)

		return self

	def run(self, cmd: str, wait=False, input=None):
		try:
			cmd = shlex.split(cmd)
			proc = PopenWithInput(cmd, shell=False, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

			if (wait):
				proc.wait()

			(result, error) = proc.communicate(input=input)

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
					os.chmod(fdir, self._mod)

				if isinstance(sfile, str):
					if(os.path.isfile(sfile)):
						if os.path.exists(os.path.join(sdir, fdir)):
							shutil.move(sfile, os.path.join(sdir, fdir, sfile))
							if files[fdir] == sfile and sfile.startswith('..'): # very hackish!
								pass
							else:
								files[fdir] = os.path.join(sdir, fdir, sfile)
						else:
							raise

				elif isinstance(sfile, list):
					for ssfile in sfile:
						if(os.path.isfile(ssfile)):
							if os.path.exists(os.path.join(sdir, fdir)):
								shutil.move(ssfile, os.path.join(sdir, fdir, ssfile))
								files[fdir][sfile.index(ssfile)] = os.path.join(sdir, fdir, ssfile)
							else:
								raise
						else:
							raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), ssfile)
				else:
					raise TypeError(f'store : [ERROR]: sfile = {type(sfile)}, must be str or list')

		except Exception:
			raise

		return files

	def clean(self):
		if not self._keep:
			if os.path.isdir(self._fdir):
				shutil.rmtree(self._fdir, ignore_errors=True)
