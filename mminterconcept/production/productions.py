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
    nsteps		= 1000	; 2 * 500000 = 1000 ps (1 ns)
    dt		    = 0.001		; 2 fs
    ; Output control
    nstxout		        = 10		; save coordinates every 10.0 ps
    nstvout		        = 10		; save velocities every 10.0 ps
    nstenergy	        = 10		; save energies every 10.0 ps
    nstlog		        = 10    	; update log file every 10.0 ps
    nstxout-compressed  = 10      ; save compressed coordinates every 10.0 ps
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
        # Yes, this isnt correct, for the illustrative mock up though, it makes it work, so meh.
        mm.LocalEnergyMinimizer.minimize(context)  # , maxIterations=500)
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

class OpenMMProduction2(OMMGeneral):

    def _make_openmm_context(self, system: mm.System) -> mm.Context:
        integrator = mm.LangevinIntegrator(298*units.kelvin,
                                           5/units.picoseconds,
                                           1 * units.femtoseconds)
        context = mm.Context(system, integrator)
        context.setPositions(self.trajectory.xyz[-1])
        context.setPeriodicBoxVectors(*self.trajectory.unitcell_vectors[-1])
        context.setVelocitiesToTemperature(298*units.kelvin)
        # Yes, this isnt correct, for the illustrative mock up though, it makes it work, so meh.
        mm.LocalEnergyMinimizer.minimize(context)  # , maxIterations=500)
        return context

    @staticmethod
    def _forces(context):
        n = context.getSystem().getNumForces()
        out = {}
        for idx in range(n):
            out[idx] = context.getState(getForces=True, groups={idx}).getForces(asNumpy=True)
        return out

    def _debug(self, context):
        state = context.getState(getPositions=True, getForces=True, getEnergy=True)
        system = context.getSystem()
        npart = system.getNumParticles()
        coords = state.getPositions(asNumpy=True)
        box = state.getPeriodicBoxVectors(asNumpy=True)
        forces = self._forces(context)
        max_forces = [f.max()._value for f in forces.values()]
        max_force = np.argmax(max_forces)
        max_particle = np.where(forces[max_force] == forces[max_force].max())[0][0]
        dist = np.zeros(npart)
        for i in range(npart):
            dist[i] = np.mean(np.sqrt(np.sum((coords[i, :] - coords[max_particle, :]) ** 2)))
        dist[max_particle] = np.inf
        return forces, coords, box, max_force, max_particle, dist

    def _make_trajectory(self, context: mm.Context) -> mdtraj.Trajectory:
        steps = 1000
        interval = 10
        frames = steps // interval
        coords = np.zeros([frames, self.trajectory.n_atoms, 3], dtype=float)
        box = np.zeros([frames, 3, 3], dtype=float)
        forces, coords_d, box_d, max_force, max_particle, dist = self._debug(context)
        for frame in range(frames):
            state = context.getState(getPositions=True)
            coords[frame, :, :] = state.getPositions(asNumpy=True)/units.nanometers
            box[frame, :, :] = state.getPeriodicBoxVectors(asNumpy=True)/units.nanometers
            print(frame)
            context.getIntegrator().step(interval)
        traj = mdtraj.Trajectory(coords, topology=self.trajectory.top)
        traj.unitcell_vectors = box
        return traj
