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

from .analysis import (analysis,
                       mda_com,
                       mda_density,
                       mda_rdf,
                       mda_rmsd,
                       mda_rog,
                       mdtraj_com,
                       mdtraj_density,
                       mdtraj_rdf,
                       mdtraj_rmsd,
                       mdtraj_rog)


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    pdbID = '1LFH'
    mdp = os.path.join(os.path.realpath(__file__), "data", "em.mdp")
    wdir = os.getcwd()

    # Construction steps
    trajectory, topology = protein_constructor(pdbID, mdp, wdir=wdir)

    # Minimization
    trajectory, topology = min_eq(trajectory=trajectory, topology=topology)

    # Equilibration
    trajectory, topology = min_eq(trajectory=trajectory, topology=topology)

    # Production
    trajectory = production(trajectory=trajectory, topology=topology)

    # Analysis
    com = analysis(trajectory)
    density = analysis(trajectory)
    rdf_r, rdf_p = analysis(trajectory)
    rmsd = analysis(trajectory)
    rog = analysis(trajectory)
