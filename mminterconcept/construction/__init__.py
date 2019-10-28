from .gmx_api import Gmx
from .example import initialize, solvate, ionize
import mdtraj
from typing import Tuple

__all__ = ["protein_constructor",
           "protein_constructor_vacuum",
           "protein_constructor_solvate_ion",
           "protein_constructor_solvate"
           ]


def protein_constructor(pdbID, box, mdp, **kwargs) -> Tuple[mdtraj.Trajectory, str]:
    raise NotImplementedError("Generic constructor does not do anything yet")


def protein_constructor_vacuum(pdbID, box, mdp, **kwargs):
    Eng = initialize(mdp, pdbID=pdbID, box=box, **kwargs)
    return Eng.getSystem(), Eng.getTop()


def protein_constructor_solvate_ion(pdbID, box, mdp, **kwargs):
    """ Returns mdtraj.Trajectory """
    if "salinity" in kwargs:
        salinity = kwargs["salinity"]
    else:
        salinity = 0.1
    Eng = initialize(mdp, pdbID=pdbID, box=box, **kwargs)
    Eng = solvate(Eng, mdp)
    Eng = ionize(Eng, salinity, mdp)

    return Eng.getSystem(), Eng.getTop()


def protein_constructor_solvate(pdbID, box, mdp, **kwargs):
    Eng = initialize(mdp, pdbID=pdbID, box=box, **kwargs)
    Eng = solvate(Eng, mdp)

    return Eng.getSystem(), Eng.getTop()