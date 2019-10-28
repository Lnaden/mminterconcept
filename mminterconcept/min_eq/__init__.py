import mdtraj
from .gro_min_eq import GroSD, GroNVT, GroNPT
from .openmm_min_eq import OpenMMMin, OpenMMEq
from typing import Tuple

__all__ = ["min_eq",
           "gro_minimize",
           "gro_eq_nvt",
           "gro_eq_npt",
           "omm_minimize",
           "omm_eq"]


def min_eq(*, trajectory: mdtraj.Trajectory, topology: str) -> Tuple[mdtraj.Trajectory, str]:
    raise NotImplementedError("Abstract model has no implementation")


def gro_minimize(*, trajectory: mdtraj.Trajectory, topology: str) -> Tuple[mdtraj.Trajectory, str]:
    return GroSD(trajectory=trajectory, topology=topology).run()


def gro_eq_nvt(*, trajectory: mdtraj.Trajectory, topology: str) -> Tuple[mdtraj.Trajectory, str]:
    return GroNVT(trajectory=trajectory, topology=topology).run()


def gro_eq_npt(*, trajectory: mdtraj.Trajectory, topology: str) -> Tuple[mdtraj.Trajectory, str]:
    return GroNPT(trajectory=trajectory, topology=topology).run()


def omm_minimize(*, trajectory: mdtraj.Trajectory, topology: str) -> Tuple[mdtraj.Trajectory, str]:
    return OpenMMMin(trajectory=trajectory, topology=topology).run()


def omm_eq(*, trajectory: mdtraj.Trajectory, topology: str) -> Tuple[mdtraj.Trajectory, str]:
    return OpenMMEq(trajectory=trajectory, topology=topology).run()

