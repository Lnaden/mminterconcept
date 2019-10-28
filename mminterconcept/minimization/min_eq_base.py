from abc import ABC, abstractmethod
import mdtraj
from pydantic import BaseModel, validator
from typing import Dict, Any

from qcelemental.models.types import Array
from qcengine.util import temporary_directory
from simtk.openmm.app.gromacstopfile import GromacsTopFile as GroTop
from contextlib import contextmanager
import os

@contextmanager
def tempcd(*args, **kwargs):
    cur = os.getcwd()
    try:
        with temporary_directory(*args, **kwargs) as tmpdir:
            os.chdir(tmpdir)
            yield
    finally:
        os.chdir(cur)


class Worklet(BaseModel, ABC):

    _output_model = BaseModel

    @abstractmethod
    def _run_internal(self) -> Dict[str, Any]:
        pass

    def run(self) -> Any:
        output = self._output_model(**self._run_internal())
        return [getattr(output, field) for field in output.__fields__.keys()]


class Trajectory(mdtraj.Trajectory):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate_traj

    @classmethod
    def validate_traj(cls, v):
        if not isinstance(v, mdtraj.Trajectory):
            raise ValueError("Not an MDTraj Topology")
        return v


class OutputSystem(BaseModel):

    @validator("topology")
    def top_is_gro(cls, v):
        with tempcd():
            gro = "gro.top"
            with open(gro, 'w') as f:
                f.write(v)
            try:
                GroTop(gro, includeDir="/Users/levinaden/miniconda3/envs/mminter/include/gromacs/")
            except:
                raise ValueError("Topology could not be processed, ensure its a valid gromacs topology string")
        return v

    trajectory: Trajectory
    topology: str


class MinimizationEquilibration(Worklet):

    @validator("topology")
    def top_is_gro(cls, v):
        with tempcd():
            gro = "gro.top"
            with open(gro, 'w') as f:
                f.write(v)
            try:
                GroTop(gro, includeDir="/Users/levinaden/miniconda3/envs/mminter/include/gromacs/")
            except:
                raise ValueError("Topology could not be processed, ensure its a valid gromacs topology string")
        return v

    parameters: Dict[str, Any] = {}
    trajectory: Trajectory
    topology: str
    _output_model = OutputSystem
