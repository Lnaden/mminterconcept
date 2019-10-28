from min_eq_base import MinimizationEquilibration, tempcd
from typing import Dict, Any
import mdtraj
import numpy as np
from textwrap import dedent

from qcengine.util import execute, temporary_directory


class GroMinEQ(MinimizationEquilibration):

    _mdp = """"""

    def _run_internal(self) -> Dict[str, Any]:
        with temporary_directory():
            # Write coordinates
            with tempcd():
                gro_name = "grocoords.gro"
                gro_file = mdtraj.formats.GroTrajectoryFile(gro_name, mode='w')
                translated_coords = self.trajectory.xyz[-1][np.newaxis, :]
                translated_cell = self.trajectory.unitcell_vectors[-1][np.newaxis, :]
                gro_file.write(translated_coords, self.trajectory.top, unitcell_vectors=translated_cell, time=[0])
                gro_file.close()
                with open(gro_name, 'r') as f:
                    gro_coords = f.read()
            topology_name = "gro.top"
            mdp_name = "mdp.mdp"
            gmxrun_base = "gmxrun"
            grompp_tpr_name = gmxrun_base + ".tpr"
            grompp_proc = ["gmx", "grompp",
                           "-f", mdp_name,
                           "-c", gro_name,
                           "-p", topology_name,
                           "-o", grompp_tpr_name,
                           "-maxwarn", "1"
                           ]
            ret_grompp, proc_grompp = execute(grompp_proc,
                                              infiles={mdp_name: self._mdp,
                                                       gro_name: gro_coords,
                                                       topology_name: self.topology,
                                                       },
                                              outfiles=[grompp_tpr_name],
                                              as_binary=[grompp_tpr_name])

            mdrun_proc = ["gmx", "mdrun",
                          "-deffnm", gmxrun_base]
            output_traj = gmxrun_base + ".trr"
            ret, proc = execute(mdrun_proc,
                                infiles={grompp_tpr_name: proc_grompp['outfiles'][grompp_tpr_name]},
                                outfiles=[output_traj],
                                as_binary=[grompp_tpr_name, output_traj])

            final_topology = self.topology
            with tempcd():
                with open(output_traj, 'wb') as f:
                    f.write(proc['outfiles'][output_traj])
                final_traj = mdtraj.load(output_traj, top=self.trajectory.topology)
            return {"trajectory": final_traj, "topology": final_topology}


class GroSD(GroMinEQ):

    _mdp = dedent("""
            ; minim.mdp - used as input into grompp to generate em.tpr
            integrator	= steep		; Algorithm (steep = steepest descent minimization)
            emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
            emstep      = 0.01      ; Energy step size
            nsteps		= 500	  	; Maximum number of (minimization) steps to perform
            
            ; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
            nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
            cutoff-scheme   = Verlet
            ns_type		    = grid		; Method to determine neighbor list (simple, grid)
            coulombtype	    = PME		; Treatment of long range electrostatic interactions
            rcoulomb	    = 1.0		; Short-range electrostatic cut-off
            rvdw		    = 1.0		; Short-range Van der Waals cut-off
            pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)
            """)


class GroNVT(GroMinEQ):

    _mdp = dedent("""
            title		= NVT equilibration 
            define		= -DPOSRES	; position restrain the protein
            ; Run parameters
            integrator	= md		; leap-frog integrator
            nsteps		= 500		; 2 * 50000 = 100 ps
            dt		    = 0.002		; 2 fs
            ; Output control
            nstxout		= 500		; save coordinates every 1.0 ps
            nstvout		= 500		; save velocities every 1.0 ps
            nstenergy	= 500		; save energies every 1.0 ps
            nstlog		= 500		; update log file every 1.0 ps
            ; Bond parameters
            continuation	        = no		; first dynamics run
            constraint_algorithm    = lincs	    ; holonomic constraints 
            constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
            lincs_iter	            = 1		    ; accuracy of LINCS
            lincs_order	            = 4		    ; also related to accuracy
            ; Neighborsearching
            cutoff-scheme   = Verlet
            ns_type		    = grid		; search neighboring grid cells
            nstlist		    = 10		; 20 fs, largely irrelevant with Verlet
            rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
            rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
            ; Electrostatics
            coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics
            pme_order	    = 4		; cubic interpolation
            fourierspacing	= 0.16	; grid spacing for FFT
            ; Temperature coupling is on
            tcoupl		= V-rescale	            ; modified Berendsen thermostat
            tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
            tau_t		= 0.1	  0.1           ; time constant, in ps
            ref_t		= 300 	  300           ; reference temperature, one for each group, in K
            ; Pressure coupling is off
            pcoupl		= no 		; no pressure coupling in NVT
            ; Periodic boundary conditions
            pbc		= xyz		    ; 3-D PBC
            ; Dispersion correction
            DispCorr	= EnerPres	; account for cut-off vdW scheme
            ; Velocity generation
            gen_vel		= yes		; assign velocities from Maxwell distribution
            gen_temp	= 300		; temperature for Maxwell distribution
            gen_seed	= -1		; generate a random seed
    """)


class GroNPT(GroMinEQ):

    _mdp = dedent("""
            title		= OPLS Lysozyme NPT equilibration 
            define		= -DPOSRES	; position restrain the protein
            ; Run parameters
            integrator	= md		; leap-frog integrator
            nsteps		= 500		; 2 * 50000 = 100 ps
            dt		    = 0.002		; 2 fs
            ; Output control
            nstxout		= 500		; save coordinates every 1.0 ps
            nstvout		= 500		; save velocities every 1.0 ps
            nstenergy	= 500		; save energies every 1.0 ps
            nstlog		= 500		; update log file every 1.0 ps
            ; Bond parameters
            continuation	        = yes		; Restarting after NVT 
            constraint_algorithm    = lincs	    ; holonomic constraints 
            constraints	            = all-bonds	; all bonds (even heavy atom-H bonds) constrained
            lincs_iter	            = 1		    ; accuracy of LINCS
            lincs_order	            = 4		    ; also related to accuracy
            ; Neighborsearching
            cutoff-scheme   = Verlet
            ns_type		    = grid		; search neighboring grid cells
            nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
            rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
            rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
            ; Electrostatics
            coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
            pme_order	    = 4		    ; cubic interpolation
            fourierspacing	= 0.16		; grid spacing for FFT
            ; Temperature coupling is on
            tcoupl		= V-rescale	            ; modified Berendsen thermostat
            tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
            tau_t		= 0.1	  0.1	        ; time constant, in ps
            ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
            ; Pressure coupling is on
            pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
            pcoupltype	        = isotropic	            ; uniform scaling of box vectors
            tau_p		        = 2.0		            ; time constant, in ps
            ref_p		        = 1.0		            ; reference pressure, in bar
            compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
            refcoord_scaling    = com
            ; Periodic boundary conditions
            pbc		= xyz		; 3-D PBC
            ; Dispersion correction
            DispCorr	= EnerPres	; account for cut-off vdW scheme
            ; Velocity generation
            gen_vel		= no		; Velocity generation is off 
    """)
