from .productions import GromacsProduction, OpenMMProduction, OpenMMProduction2
from typing import Tuple
import mdtraj

__all__ = [
    "production",
    "gromacs_production",
    "openmm_production",
    "openmm_production2"
]


def production(*, trajectory: mdtraj.Trajectory, topology: str) -> mdtraj.Trajectory:
    raise NotImplementedError("Abstract model has no implementation")


def gromacs_production(*, trajectory: mdtraj.Trajectory, topology: str) -> mdtraj.Trajectory:
    traj, _ = GromacsProduction(trajectory=trajectory, topology=topology).run()
    return traj


def openmm_production(*, trajectory: mdtraj.Trajectory, topology: str) -> mdtraj.Trajectory:
    traj, _ = OpenMMProduction(trajectory=trajectory, topology=topology).run()
    return traj

def openmm_production2(*, trajectory: mdtraj.Trajectory, topology: str) -> mdtraj.Trajectory:
    traj, _ = OpenMMProduction2(trajectory=trajectory, topology=topology).run()
    return traj
