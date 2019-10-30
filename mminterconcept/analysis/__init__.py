from .models import Component
from .mda_com import COMMDAnalysisComponent
from .mda_density import DensityMDAnalysisComponent
from .mda_rdf import RDFMDAnalysisComponent
from .mda_rmsd import RMSDMDAnalysisComponent
from .mda_rog import ROGMDAnalysisComponent

from .mdtraj_com import COMComponent
from .mdtraj_density import DensityComponent
from .mdtraj_rdf import RDFComponent
from .mdtraj_rmsd import RMSDComponent
from .mdtraj_rog import ROGComponent

import mdtraj
import numpy as np

from typing import Tuple, Union
from .plot import *

__all__ = ["analysis",
           "mda_com",
           "mda_den",
           "mda_rdf",
           "mda_rmsd",
           "mda_rog",
           "mdtraj_com",
           "mdtraj_den",
           "mdtraj_rdf",
           "mdtraj_rmsd",
           "mdtraj_rog"
           ]


def analysis(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    raise NotImplementedError()

# COM
def mda_com(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = COMMDAnalysisComponent(top, trajectory)
    return a.run(top, trajectory)

def mdtraj_com(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = COMComponent(top, trajectory)
    return a.run(top, trajectory)

# Density
def mda_den(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = DensityMDAnalysisComponent(top, trajectory)
    return a.run(top, trajectory)

def mdtraj_den(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = DensityComponent(top, trajectory)
    return a.run(top, trajectory)

# RDF
def mda_rdf(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFMDAnalysisComponent(top, trajectory)
    return a.run(top, trajectory)

def mdtraj_rdf(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFComponent(top, trajectory)
    return a.run(top, trajectory)

# RMSD
def mda_rmsd(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = RMSDMDAnalysisComponent(top, trajectory)
    return a.run(top, trajectory)

def mdtraj_rmsd(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = RMSDComponent(top, trajectory)
    return a.run(top, trajectory)

# RoG
def mda_rog(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = ROGMDAnalysisComponent(top, trajectory)
    return a.run(top, trajectory)

def mdtraj_rog(top: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = ROGComponent(top, trajectory)
    return a.run(top, trajectory)
