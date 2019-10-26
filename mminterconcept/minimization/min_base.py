from abc import ABC, abstractmethod
import mdtraj
from pydantic import BaseModel
from typing import Dict, Any, Union

from qcelemental.models.types import Array


class Worklet(BaseModel, ABC):

    _output_model = BaseModel

    @abstractmethod
    def _run_internal(self) -> Dict[str, Any]:
        pass

    def run(self) -> Any:
        return self._output_model(**self._run_internal())


class Component:
    def __init__(self, input_model, output_model):
        self.input_model = input_model
        self.output_model = output_model

    def __call__(self, **kwargs):
        model_init = self.input_model(**kwargs)


class Topology(mdtraj.Topology):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate_top

    @classmethod
    def validate_top(cls, v):
        if not isinstance(v, mdtraj.Topology):
            raise ValueError("Not an MDTraj Topology")
        return v


class System(BaseModel):
    topology: Union[Topology, str]
    coordinates: Array[float]
    unit_cell: Array[float] = None


class ForceField(BaseModel):
    everything: str = None
    solute: str = None
    solvent: str = None


class Minimization(Worklet):
    parameters: Dict[str, Any] = {}
    system: System
    forcefield: ForceField
    _output_model = System
