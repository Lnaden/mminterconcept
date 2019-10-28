import mdtraj
from gro_min_eq import GroSD, GroNVT, GroNPT
from openmm_min_eq import OpenMMMin, OpenMMEq

water = "tip3p"
ff = "amber99"

coords = "../data/1LFH/gmx/struct/system_ionized.gro"
top = "../data/1LFH/gmx/top/system_ionized.top"
inc = "/Users/levinaden/miniconda3/envs/mminter/include/gromacs/"

minimizer = GroSD
equlibrators = [GroNVT, GroNPT]

minimizer = OpenMMMin
# equlibrators = [GroNVT, GroNPT, OpenMMEq]
equlibrators = [OpenMMEq]


def minimize(*, trajectory, topology):
    return minimizer(trajectory=trajectory, topology=topology).run()


def equilibrate(*, trajectory, topology):
    for eq in equlibrators:
        print(f"Running Equilibration {eq}")
        trajectory, topology = eq(trajectory=trajectory, topology=topology).run()
    return trajectory, topology


with open("solvated.top", 'r') as f:
    top = f.read()

min_traj, min_top = minimize(trajectory=mdtraj.load("ion.h5"), topology=top)
eq_traj, eq_top = equilibrate(trajectory=min_traj, topology=min_top)
