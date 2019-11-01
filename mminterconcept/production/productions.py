from textwrap import dedent

import simtk.openmm as mm
import simtk.unit as units
from ..min_eq.gro_min_eq import GroMinEQ
from ..min_eq.openmm_min_eq import OMMGeneral

import mdtraj
import numpy as np


class GromacsProduction(GroMinEQ):
    _mdp = dedent("""
    title		= OPLS Lysozyme MD simulation 
    ; Run parameters
    integrator	= md		; leap-frog integrator
    nsteps		= 5000	; 2 * 500000 = 1000 ps (1 ns)
    dt		    = 0.002		; 2 fs
    ; Output control
    nstxout		        = 5000		; save coordinates every 10.0 ps
    nstvout		        = 5000		; save velocities every 10.0 ps
    nstenergy	        = 5000		; save energies every 10.0 ps
    nstlog		        = 5000		; update log file every 10.0 ps
    nstxout-compressed  = 5000      ; save compressed coordinates every 10.0 ps
                                    ; nstxout-compressed replaces nstxtcout
    compressed-x-grps   = System    ; replaces xtc-grps
    ; Bond parameters
    continuation	        = yes		; Restarting after NPT 
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
    ; Periodic boundary conditions
    pbc		= xyz		; 3-D PBC
    ; Dispersion correction
    DispCorr	= EnerPres	; account for cut-off vdW scheme
    ; Velocity generation
    gen_vel		= no		; Velocity generation is off 
    """)


class OpenMMProduction(OMMGeneral):

    def _make_openmm_context(self, system: mm.System) -> mm.Context:
        integrator = mm.LangevinIntegrator(298*units.kelvin,
                                           5/units.picoseconds,
                                           1 * units.femtoseconds)
        context = mm.Context(system, integrator)
        context.setPositions(self.trajectory.xyz[-1])
        context.setPeriodicBoxVectors(*self.trajectory.unitcell_vectors[-1])
        context.setVelocitiesToTemperature(298*units.kelvin)
        return context

    def _make_trajectory(self, context: mm.Context) -> mdtraj.Trajectory:
        steps = 1000
        interval = 10
        frames = steps // interval
        coords = np.zeros([frames, self.trajectory.n_atoms, 3], dtype=float)
        box = np.zeros([frames, 3, 3], dtype=float)
        for frame in range(frames):
            state = context.getState(getPositions=True)
            coords[frame, :, :] = state.getPositions(asNumpy=True)/units.nanometers
            box[frame, :, :] = state.getPeriodicBoxVectors(asNumpy=True)/units.nanometers
            context.getIntegrator().step(interval)
        traj = mdtraj.Trajectory(coords, topology=self.trajectory.top)
        traj.unitcell_vectors = box
        return traj
