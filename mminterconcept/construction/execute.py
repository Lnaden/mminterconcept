'''
Author: Andrew-AbiMansour
DOC: Oct 22, 2019

Gromacs API
'''

import subprocess
from tools import RandString
import os

class Process:
	def __init__(self, files=None):
		self.files = files

	def __enter__(self):
		self.fdir = RandString.name()
		if os.path.isdir(self.fdir):
			raise Exception(f'__enter__ : [ERROR]: randrom dir {self.fdir} already exists')
		else:
			os.mkdir(self.fdir)
			os.chdir(self.fdir)
			print(os.getcwd())

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

	def clean(self):
		if isinstance(self.files, list):
			for file in self.files:
				if(os.path.exists(file)):
					os.remove(file)
		elif self.files:
			if os.path.exists(self.files):
				os.remove(self.files)

		if os.path.isdir(self.fdir):
			os.rmdir(self.fdir)

with Process() as fp:
	fp.run(cmd=['ls', '.'])
