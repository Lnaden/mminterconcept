"""
mminterconcept.py
The concept implementation for the MolSSI MM Interoperable Components Project

Handles the primary functions
"""

import os
from .construction import (protein_constructor,
                           protein_constructor_vacuum,
                           protein_constructor_solvate,
                           protein_constructor_solvate_ion)

from .min_eq import (min_eq,
                     gro_minimize,
                     gro_eq_nvt,
                     gro_eq_npt,
                     omm_minimize,
                     omm_eq)

from .production import (production,
                         gromacs_production,
                         openmm_production)


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    pdbID = '1LFH'
    box = (9.5, 9, 7)
    mdp = os.path.join(os.path.realpath(__file__), "data", "em.mdp")
    wdir = os.getcwd()

    # Construction steps
    trajectory, topology = protein_constructor(pdbID, box, mdp, wdir=wdir)

    # Minimization
    trajectory, topology = min_eq(trajectory=trajectory, topology=topology)

    # Equilibration
    trajectory, topology = min_eq(trajectory=trajectory, topology=topology)

    # Production
    trajectory = production(trajectory=trajectory, topology=topology)

    # Analysis

