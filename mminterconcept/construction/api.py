'''
Author: Andrew-AbiMansour
DOC: Oct 22, 2019

Generic API for MD solvers (e.g. Gromacs)
'''

import subprocess
from tools import RandString, ImportPDB
from core import Box
from collections import OrderedDict
import os, sys
import errno
import shutil, shlex
import six

class PopenWithInput(subprocess.Popen):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def communicate(self, input=None):
		if input:
			return super().communicate(input)
		else:
			return super().communicate()

class Engine:
	def __init__(self, pdbID : str, ff_solute: str, ff_solvent: str, topfname: str, ofname: str, ext: str, exec: str, **args):

		if 'fdir' in args:
			self._fdir = args['fdir']
		else:
			self._fdir = RandString.name()

		if 'mod' not in args:
			self._mod = 0o777

		if os.path.isdir(self._fdir):
			raise Exception(f'__enter__ : [ERROR]: randrom dir {self._fdir} already exists')
		else:
			os.mkdir(self._fdir)
			os.chmod(self._fdir, self._mod)

		os.chdir(self._fdir)

		self._args = {
			'_exec': exec,
			'ifname': os.path.abspath(f'{pdbID}.{ext}'),
			'ff_solute': ff_solute,
			'ff_solvent': ff_solvent,
			'ofname': ofname,
			'topfname': topfname
		}

		if 'sdir' in args:
			self._sdir = sdir
		else:
			self._sdir = '..'

		if 'box' not in args:
			self._Box = Box(shape='cubic', bound=(10.0, 10.0, 10.0)) # need to fill this from pdb file
		else:
			self._Box = Box(shape='cubic', bound=args['box'])

		self._pdbID = pdbID

		# dict not always ordered (depends on python ver)
		self._gmx_args = OrderedDict()

		with ImportPDB(self._pdbID) as SS:
			SS.save(f'{pdbID}.{ext}')

		if not os.path.exists('struct'):
			os.mkdir('struct')
			os.chmod('struct', self._mod)

		if not os.path.exists('top'):
			os.mkdir('top')
			os.chmod('top', self._mod)

		if not os.path.exists('tpr'):
			os.mkdir('tpr')
			os.chmod('tpr', self._mod)

	def _dict_to_str(self, **args):
		""" Map dict to flattened list """
		args = [[key, val] for key, val in args.items()]
		args = [single for pair in args for single in pair]
		return ' '.join([arg for arg in args if not arg.startswith('_')])

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
			proc = PopenWithInput(cmd)

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
							print(sdir, fdir)
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

	def Popen(self, *args, **kwargs):
		"""Returns a special Popen instance (:class:`PopenWithInput`). """

		stderr = kwargs.pop('stderr', None)     # default: print to stderr (if STDOUT then merge)

		if stderr is False:                     # False: capture it
			stderr = PIPE
		elif stderr is True:
			stderr = None                       # use stderr

		stdout = kwargs.pop('stdout', None)     # either set to PIPE for capturing output

		if stdout is False:                     # ... or to False
			stdout = PIPE
		elif stdout is True:
			stdout = None                       # for consistency, make True write to screen

		stdin = kwargs.pop('stdin', None)
		input = kwargs.pop('input', None)

		use_shell = kwargs.pop('use_shell', False)

		if input:
			stdin = PIPE
		if isinstance(input, six.string_types) and not input.endswith('\n'):
			# make sure that input is a simple string with \n line endings
			input = six.text_type(input) + '\n'
		else:
			try:
				# make sure that input is a simple string with \n line endings
				input = '\n'.join(map(six.text_type, input)) + '\n'
			except TypeError:
				# so maybe we are a file or something ... and hope for the best
				pass

		cmd = self._commandline(*args, **kwargs)   # lots of magic happening here
				# (cannot move out of method because filtering of stdin etc)
		try:
			proc = PopenWithInput(cmd, stdin=stdin, stderr=stderr, stdout=stdout,
				universal_newlines=True, input=input, shell=use_shell)
		except: raise

		if err.errno == errno.ENOENT:
			errmsg = "Failed to find Gromacs command {0!r}, maybe its not on PATH or GMXRC must be sourced?".format(self.command_name)
			raise OSError(errmsg)
		else:
			raise

		return proc
