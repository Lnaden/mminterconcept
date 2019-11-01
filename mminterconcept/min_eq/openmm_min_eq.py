from .min_eq_base import MinimizationEquilibration, tempcd
from typing import Dict, Any
import mdtraj
import parmed
import numpy as np
import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as units

from abc import abstractmethod

GroTop = app.gromacstopfile.GromacsTopFile


class OMMMinEq(MinimizationEquilibration):

    @abstractmethod
    def _do_openmm_thing(self, system: mm.System) -> mm.Context:
        pass

    def _run_internal(self) -> Dict[str, Any]:
        # Create System
        with tempcd():
            gro = "gro.top"
            with open(gro, 'w') as f:
                f.write(self.topology)
            # Convert gromacs top to parmed top before reading BECAUSE REASONS?
            p_gro = parmed.gromacs.GromacsTopologyFile(gro,
                                                       # Box is in angstroms
                                                       box=[*(self.trajectory.unitcell_lengths * 10)[-1],
                                                            *self.trajectory.unitcell_angles[-1]],
                                                       parametrize=True)

            # g_top = GroTop(gro,
            #                periodicBoxVectors=self.trajectory.unitcell_vectors[-1],
            #                includeDir="/Users/levinaden/miniconda3/envs/mminter/include/gromacs/")
            # system = g_top.createSystem(nonbondedMethod=app.forcefield.PME,
            #                             #constraints=app.forcefield.HBonds,
            #                             constraints=app.forcefield.AllBonds,
            #                             rigidWater=True,
            #                             nonbondedCutoff=1 * units.nanometer)
            system = p_gro.createSystem(nonbondedMethod=app.forcefield.PME,
                                        nonbondedCutoff=1 * units.nanometer,
                                        constraints=app.forcefield.AllBonds,
                                        rigidWater=True
                                        )
            # system.addForce(mm.AndersenThermostat(298 * units.kelvin, 1 / units.picosecond))
            has_barostat = False
            for force in system.getForces():
                if isinstance(force, mm.MonteCarloBarostat):
                    has_barostat = True
                    break
            if not has_barostat:
                system.addForce(mm.MonteCarloBarostat(1*units.bar, 298*units.kelvin))
            output_context = self._do_openmm_thing(system)
            final_state = output_context.getState(getPositions=True, enforcePeriodicBox=False)
            final_coords = final_state.getPositions(asNumpy=True)
            final_vector = final_state.getPeriodicBoxVectors(asNumpy=True)
            output_traj = mdtraj.Trajectory(final_coords[np.newaxis, :]/units.nanometers, topology=self.trajectory.top)
            output_traj.unitcell_vectors = final_vector[np.newaxis, :]/units.nanometers
        return {"trajectory": output_traj, "topology": self.topology}


# Wholesale borrowed from OpenMMTools and John Chodera
class FIREMinimizationIntegrator(mm.CustomIntegrator):
    """Fast Internal Relaxation Engine (FIRE) minimization.
    Notes
    -----
    This integrator is taken verbatim from Peter Eastman's example appearing in the CustomIntegrator header file documentation.
    References
    ----------
    Erik Bitzek, Pekka Koskinen, Franz Gaehler, Michael Moseler, and Peter Gumbsch.
    Structural Relaxation Made Simple. PRL 97:170201, 2006.
    http://dx.doi.org/10.1103/PhysRevLett.97.170201
    Examples
    --------
    Create a FIRE integrator with default parameters.
    >>> integrator = FIREMinimizationIntegrator()
    """

    def __init__(self, timestep=1.0 * units.femtoseconds, tolerance=None, alpha=0.1, dt_max=10.0 * units.femtoseconds, f_inc=1.1, f_dec=0.5, f_alpha=0.99, N_min=5):
        """Construct a Fast Internal Relaxation Engine (FIRE) minimization integrator.
        Parameters
        ----------
        timestep : simtk.unit.Quantity compatible with femtoseconds, optional, default = 1*femtoseconds
            The integration timestep.
        tolerance : simtk.unit.Quantity compatible with kilojoules_per_mole/nanometer, optional, default = None
            Minimization will be terminated when RMS force reaches this tolerance.
        alpha : float, optional default = 0.1
            Velocity relaxation parameter, alpha \in (0,1).
        dt_max : simtk.unit.Quantity compatible with femtoseconds, optional, default = 10*femtoseconds
            Maximum allowed timestep.
        f_inc : float, optional, default = 1.1
            Timestep increment multiplicative factor.
        f_dec : float, optional, default = 0.5
            Timestep decrement multiplicative factor.
        f_alpha : float, optional, default = 0.99
            alpha multiplicative relaxation parameter
        N_min : int, optional, default = 5
            Limit on number of timesteps P is negative before decrementing timestep.
        Notes
        -----
        Velocities should be set to zero before using this integrator.
        """

        # Check input ranges.
        if not ((alpha > 0.0) and (alpha < 1.0)):
            raise Exception("alpha must be in the interval (0,1); specified alpha = %f" % alpha)

        if tolerance is None:
            tolerance = 0 * units.kilojoules_per_mole / units.nanometers

        super(FIREMinimizationIntegrator, self).__init__(timestep)

        self.addGlobalVariable("alpha", alpha)  # alpha
        self.addGlobalVariable("P", 0)  # P
        self.addGlobalVariable("N_neg", 0.0)
        self.addGlobalVariable("fmag", 0)  # |f|
        self.addGlobalVariable("fmax", 0)  # max|f_i|
        self.addGlobalVariable("ndof", 0)  # number of degrees of freedom
        self.addGlobalVariable("ftol", tolerance.value_in_unit_system(units.md_unit_system))  # convergence tolerance
        self.addGlobalVariable("freg", 0.001)  # regularization to add to |f|
        self.addGlobalVariable("vmag", 0)  # |v|
        self.addGlobalVariable("converged", 0) # 1 if convergence threshold reached, 0 otherwise
        self.addPerDofVariable("x0", 0)
        self.addPerDofVariable("x1", 0)

        # Enclose everything in a block that checks if we have already converged.
        self.beginIfBlock('converged < 1')

        # Update Context.
        # TODO: Remove this?
        self.addUpdateContextState()

        # Store old positions
        self.addComputePerDof('x0', 'x')

        # MD: Take a velocity Verlet step.
        self.addComputePerDof("v", "v+0.5*dt*f/m")
        self.addComputePerDof("x", "x+dt*v")
        self.addComputePerDof("x1", "x")
        self.addConstrainPositions()
        self.addComputePerDof("v", "v+0.5*dt*f/m+(x-x1)/dt")
        self.addConstrainVelocities()

        # Compute fmag = |f|
        self.addComputeGlobal('fmag', '0.0')
        self.addComputeSum('fmag', 'f*f')
        self.addComputeGlobal('fmag', 'sqrt(fmag)')
        # Compute vmag = |v|
        self.addComputeGlobal('vmag', '0.0')
        self.addComputeSum('vmag', 'v*v')
        self.addComputeGlobal('vmag', 'sqrt(vmag)')

        # Check for convergence of force magnitude.
        # TODO: Check componentwise or RMS?
        self.addComputeSum('ndof', '1')
        #self.addComputeSum('converged', 'step(ftol - abs(f)) / ndof')
        self.beginIfBlock('fmag/sqrt(ndof) < ftol')
        self.addComputeGlobal('converged', '1') # signal convergence
        self.endBlock()

        ## Check for NaN and try to recover.
        #self.beginIfBlock('abs(step(energy)-step(1-energy)) = 0')
        ## Reset positions and set velocities to zero.
        #self.addComputePerDof('x', 'x0')
        #self.addComputePerDof('v', '0.0')
        ## Reset P counter.
        #self.addComputeGlobal('N_neg', '0.0')
        ## Scale down timestep.
        #self.addComputeGlobal('dt', 'dt*%f' % f_dec)
        #self.endBlock()

        # F1: Compute P = F.v
        self.addComputeSum('P', 'f*v')

        # F2: set v = (1-alpha) v + alpha \hat{F}.|v|
        # Update velocities.
        self.addComputePerDof('v', '(1-alpha)*v + alpha*(f/(fmag+freg))*vmag')

        # F3: If P > 0 and the number of steps since P was negative > N_min,
        # Increase timestep dt = min(dt*f_inc, dt_max) and decrease alpha = alpha*f_alpha
        self.beginIfBlock('P > 0')
        # Update count of number of steps since P was negative.
        self.addComputeGlobal('N_neg', 'step(P) * (N_neg + 1)')
        # If we have enough steps since P was negative, scale up timestep.
        self.beginIfBlock('N_neg > %d' % N_min)
        self.addComputeGlobal('dt', 'min(dt*%f, %f)' % (f_inc, dt_max.value_in_unit_system(units.md_unit_system))) # TODO: Automatically convert dt_max to md units
        self.addComputeGlobal('alpha', 'alpha * %f' % f_alpha)
        self.endBlock()
        self.endBlock()

        # F4: If P < 0, decrease the timestep dt = dt*f_dec, freeze the system v=0,
        # and set alpha = alpha_start
        self.beginIfBlock('P < 0')
        self.addComputeGlobal('N_neg', '0.0')
        self.addComputeGlobal('dt', 'dt*%f' % f_dec)
        self.addComputePerDof('v', '0.0')
        self.endBlock()

        # Close block that checks for convergence.
        self.endBlock()


class OpenMMMin(OMMMinEq):

    @staticmethod
    def _build_integrator() -> mm.Integrator:
        # return mm.VerletIntegrator(0.1 * units.femtoseconds)
        return mm.LangevinIntegrator(298*units.kelvin,
                                     5/units.picoseconds,
                                     1 * units.femtoseconds)

    def _do_openmm_thing(self, system: mm.System) -> mm.Context:
        context = mm.Context(system, self._build_integrator(), mm.Platform.getPlatformByName("CPU"))
        context.setPositions(self.trajectory.xyz[-1])
        context.setPeriodicBoxVectors(*self.trajectory.unitcell_vectors[-1])
        mm.LocalEnergyMinimizer.minimize(context, tolerance=2)  # , maxIterations=500)
        # context.getIntegrator().step(50)
        # pos = context.getState(getPositions=True).getPositions(asNumpy=True)
        # mm.LocalEnergyMinimizer.minimize(context, tolerance=0.5)  # , maxIterations=500)
        return context


class OpenMMFireMin(OpenMMMin):

    @staticmethod
    def _build_integrator() -> mm.Integrator:
        return FIREMinimizationIntegrator()


class OpenMMEq(OMMMinEq):

    def _do_openmm_thing(self, system: mm.System) -> mm.Context:
        integrator = mm.LangevinIntegrator(298*units.kelvin,
                                           5/units.picoseconds,
                                           1 * units.femtoseconds)
        # integrator = mm.VerletIntegrator(1 * units.femtoseconds)
        context = mm.Context(system, integrator)
        context.setPositions(self.trajectory.xyz[-1])
        context.setPeriodicBoxVectors(*self.trajectory.unitcell_vectors[-1])
        context.setVelocitiesToTemperature(298*units.kelvin)
        integrator.step(10)
        return context
