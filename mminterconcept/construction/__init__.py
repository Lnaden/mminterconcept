from .gmx_api import Gmx
from .example import initialize, solvate, ionize
import mdtraj
from typing import Tuple

__all__ = ["protein_constructor",
           "protein_constructor_vacuum",
           "protein_constructor_solvate_ion",
           "protein_constructor_solvate"
           ]

def protein_constructor(pdbID, mdp, **kwargs) -> Tuple[mdtraj.Trajectory, str]:
    raise NotImplementedError("Generic constructor does not do anything yet")


def protein_constructor_vacuum(pdbID, mdp, **kwargs):
    Eng = initialize(mdp, pdbID=pdbID, **kwargs)
    Sys, top = Eng.getSystem(), Eng.getTop()
    Eng.clean()
    return Sys, top


def protein_constructor_solvate_ion(pdbID, mdp, **kwargs):
    """ Returns mdtraj.Trajectory """
    if "salinity" in kwargs:
        salinity = kwargs["salinity"]
    else:
        salinity = 0.1
    Eng = initialize(mdp, pdbID=pdbID, **kwargs)
    Eng = solvate(Eng, mdp)
    Eng = ionize(Eng, salinity, mdp)
    Sys, top = Eng.getSystem(), Eng.getTop()
    Eng.clean()
    return Sys, top

def protein_constructor_solvate(pdbID, mdp, **kwargs):
    Eng = initialize(mdp, pdbID=pdbID, **kwargs)
    Eng = solvate(Eng, mdp)
    Sys, top = Eng.getSystem(), Eng.getTop()
    Eng.clean()
    return Sys, top
