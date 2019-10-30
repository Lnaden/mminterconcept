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


__all__ = ["analysis",
           "mda_com",
           "mda_density",
           "mda_rdf",
           "mda_rmsd",
           "mda_rog",
           "mdtraj_com",
           "mdtraj_density",
           "mdtraj_rdf",
           "mdtraj_rmsd",
           "mdtraj_rog"
           ]


def analysis(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    raise NotImplementedError()

# COM
def mda_com(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = COMMDAnalysisComponent(struct, trajectory)
    return a.run()


def mdtraj_com(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = COMComponent(struct, trajectory)
    return a.run()

# Density
def mda_density(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = DensityMDAnalysisComponent(struct, trajectory)
    return a.run()

def mdtraj_density(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = DensityComponent(struct, trajectory)
    return a.run()

# RDF
def mda_rdf(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFMDAnalysisComponent(struct, trajectory)
    return a.run()

def mdtraj_rdf(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFComponent(struct, trajectory)
    return a.run()

# RMSD
def mda_rmsd(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = RMSDMDAnalysisComponent(struct, trajectory)
    return a.run()

def mdtraj_rmsd(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = RMSDComponent(struct, trajectory)
    return a.run()

# RoG
def mda_rog(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = ROGMDAnalysisComponent(struct, trajectory)
    return a.run()

def mdtraj_rog(struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = ROGComponent(struct, trajectory)
    return a.run()
