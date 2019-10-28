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
    return a.compute()


def mdtraj_com(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = COMComponent(trajectory)
    return a.compute()


# Density
def mda_density(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = DensityMDAnalysisComponent(trajectory)
    return a.compute()


def mdtraj_density(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = DensityComponent(trajectory)
    return a.compute()


# RDF
def mda_rdf(trajectory: mdtraj.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFMDAnalysisComponent(trajectory)
    return a.compute()


def mdtraj_rdf(trajectory: mdtraj.Trajectory) -> Tuple[np.ndarray, np.ndarray]:
    a = RDFComponent(trajectory)
    return a.compute()


# RMSD
def mda_rmsd(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = RMSDMDAnalysisComponent(trajectory)
    return a.compute()


def mdtraj_rmsd(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = RMSDComponent(trajectory)
    return a.compute()


# RoG
def mda_rog(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = ROGMDAnalysisComponent(trajectory)
    return a.compute()


def mdtraj_rog(trajectory: mdtraj.Trajectory) -> np.ndarray:
    a = ROGComponent(trajectory)
    return a.compute()
