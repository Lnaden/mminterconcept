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


def analysis(trajectory: mdtraj.Trajectory) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    raise NotImplementedError()


# COM
def mda_com(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = COMMDAnalysisComponent(trajectory)
    return a.run()


def mdtraj_com(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = COMComponent(trajectory)
    return a.run()


# Density
def mda_density(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = DensityMDAnalysisComponent(trajectory)
    return a.run()


def mdtraj_density(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = DensityComponent(trajectory)
    return a.run()


# RDF
def mda_rdf(trajectory: mdtraj.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFMDAnalysisComponent(trajectory)
    return a.run()


def mdtraj_rdf(trajectory: mdtraj.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFComponent(trajectory)
    return a.run()


# RMSD
def mda_rmsd(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = RMSDMDAnalysisComponent(trajectory)
    return a.run()


def mdtraj_rmsd(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = RMSDComponent(trajectory)
    return a.run()


# RoG
def mda_rog(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = ROGMDAnalysisComponent(trajectory)
    return a.run()


def mdtraj_rog(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = ROGComponent(trajectory)
    return a.run()
