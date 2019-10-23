'''
Author: Andrew-AbiMansour
DOC: Oct 22, 2019

Generic API for MD solvers (e.g. Gromacs)
'''

import subprocess
from tools import RandString
import os
import errno
import shutil

class Process:
	def __init__(self, keep=False):
		self._keep = keep

	def __enter__(self):
		self.fdir = RandString.name()
		if os.path.isdir(self.fdir):
			raise Exception(f'__enter__ : [ERROR]: randrom dir {self.fdir} already exists')
		else:
			os.mkdir(self.fdir)
			os.chdir(self.fdir)

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

	def store(self, **files):
		try:
			for fdir, file in files.items():
				if not os.path.exists(fdir):
					os.mkdir(fdir)

				if isinstance(file, str):
					if(os.path.isfile(file)):
						shutil.move(file, os.path.join(fdir, file))
				elif isinstance(file, list):
					for sfile in file:
						if(os.path.isfile(sfile)):
							shutil.move(sfile, os.path.join(fdir, sfile))
						else:
							raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)
				else:
					raise TypeError(f'store : [ERROR]: file = {type(file)}, must be str or list')
		except Exception:
			raise

	def clean(self):
		if not self._keep:
			if os.path.isdir(self.fdir):
				os.rmdir(self.fdir)
