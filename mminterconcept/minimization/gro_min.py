from min_base import Minimization
from typing import Dict, Any
import mdtraj
import numpy as np
from textwrap import dedent

from qcengine.util import execute, temporary_directory

from contextlib import contextmanager
import os

@contextmanager
def tempcd(*args, **kwargs):
    cur = os.getcwd()
    try:
        with temporary_directory(*args, **kwargs) as tmpdir:
            os.chdir(tmpdir)
            yield
    finally:
        os.chdir(cur)


with tempcd():
    # do stuff
    pass


class GroMinEQ(Minimization):

    _mdp = """"""

    def _run_internal(self) -> Dict[str, Any]:
        with temporary_directory():
            if isinstance(self.system.topology, str):
                gro_name = self.system.topology
            else:
                # Write coordinates
                gro_name = "grocoords.gro"
                gro_file = mdtraj.formats.GroTrajectoryFile(gro_name, mode='w')
                translated_coords = self.system.coordinates[np.newaxis, :]
                translated_cell = self.system.unit_cell[np.newaxis, :]
                gro_file.write(translated_coords, self.system.topology, unitcell_vectors=translated_cell, time=[0])
                gro_file.close()
            topology_name = "gro.top"
            processed_name = "groproc.gro"
            top_proc = ["gmx", "pdb2gmx",
                        "-f", gro_name,
                        "-ff", self.forcefield.solute,
                        "-water", self.forcefield.solvent,
                        "-ignh",
                        "-p", topology_name,
                        "-o", processed_name]
            with open(gro_name, 'r') as f:
                gro_data = f.read()
            breakpoint()
            ret_pdb, proc_pdb = execute(top_proc,
                                        infiles={gro_name: gro_data},
                                        outfiles=[topology_name])
            mdp_name = "mdp.mdp"
            gmxrun_base = "gmxrun"
            grompp_tpr_name = gmxrun_base + ".tpr"
            grompp_proc = ["gmx", "grompp",
                           "-f", mdp_name,
                           "-c", gro_name,
                           "-p", topology_name,
                           "-o", grompp_tpr_name
                           ]

            ret_grompp, proc_grompp = execute(grompp_proc,
                                              infiles={mdp_name: self._mdp,
                                                       gro_name: gro_data,
                                                       topology_name: proc_pdb['outfiles'][topology_name],
                                                       },
                                              outfiles=[grompp_tpr_name],
                                              as_binary=[grompp_tpr_name])

            mdrun_proc = ["gmx", "mdrun",
                          "-deffnm", gmxrun_base]
            output_traj = gmxrun_base + ".trr"
            ret, proc = execute(mdrun_proc,
                                infiles={mdp_name: self._mdp,
                                         gro_name: gro_data,
                                         topology_name: proc_grompp['outfiles'][topology_name],
                                         },
                                outfiles=[output_traj],
                                as_binary=[output_traj])

            final_topology = self.system.topology
            final_traj = mdtraj.load(output_traj, top=self.system.topology)
            final_coords = final_traj.xyz[-1]
            final_unitcell = final_traj.unitcell_vectors
            return {"topology": final_topology, "coordinates": final_coords, "unit_cell": final_unitcell}


class GroSD(GroMinEQ):

    _mdp = dedent("""
            ; minim.mdp - used as input into grompp to generate em.tpr
            integrator	= steep		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
            emstep      = 0.01      ; Energy step size
            nsteps		= 50	  	; Maximum number of (minimization) steps to perform
            
            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; Short-range electrostatic cut-off
            rvdw		    = 1.0		; Short-range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)
            """)

