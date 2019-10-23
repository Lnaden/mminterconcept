'''
Author: Andrew-AbiMansour
DOC: Oct 21, 2019

This module defines the System class for MD simulation/analysis
'''

from tools import ImportPDB
from pydantic import BaseModel, Field
import mdtraj
from typing import Any

class Species(mdtraj.Trajectory):
	@classmethod
	def __get_validators__(cls):
		yield cls.validate

	@classmethod
	def validate(cls, v):
		if not isinstance(v, str):
			raise ValueError(f'strict string: str expected not {type(v)}')
		return v

class Box(BaseModel):
	shape: str = 'box'
	bound: tuple

class MDparams(BaseModel):
	temp: float = None
	press: float = None

	def __init__(self, temp, press, **args):
		super().__init__(**args)

class EMparams(BaseModel):
        method: str = None
        tol: float = None

        def __init__(self, method, tol, **args):
                super().__init__(**args)

class Props(BaseModel):
	pH: float = None

	def __init__(self, pH, **args):
		super().__init__(**args)

class System:
	Box: 'SimBox'
	Species: 'Solute'

	def __init__(self, Box, Species, **args):
		pass
