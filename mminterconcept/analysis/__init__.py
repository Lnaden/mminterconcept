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
from .visualize import *

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


def analysis(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    raise NotImplementedError()

# COM
def mda_com(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> np.ndarray:
    a = COMMDAnalysisComponent(trajectory, top)
    return a.run(trajectory, top)

def mdtraj_com(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> np.ndarray:
    a = COMComponent(trajectory, top)
    return a.run(trajectory, top)

# Density
def mda_density(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> np.ndarray:
    a = DensityMDAnalysisComponent(trajectory, top)
    return a.run(trajectory, top)

def mdtraj_density(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> np.ndarray:
    a = DensityComponent(trajectory, top)
    return a.run(trajectory, top)

# RDF
def mda_rdf(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFMDAnalysisComponent(trajectory, top)
    return a.run(trajectory, top)

def mdtraj_rdf(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFComponent(trajectory, top)
    return a.run(trajectory, top)

# RMSD
def mda_rmsd(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> np.ndarray:
    a = RMSDMDAnalysisComponent(trajectory, top)
    return a.run(trajectory, top)

def mdtraj_rmsd(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> np.ndarray:
    a = RMSDComponent(trajectory, top)
    return a.run(trajectory, top)

# RoG
def mda_rog(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> np.ndarray:
    a = ROGMDAnalysisComponent(trajectory, top)
    return a.run(trajectory, top)

def mdtraj_rog(trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory=None) -> np.ndarray:
    a = ROGComponent(trajectory, top)
    return a.run(trajectory, top)
