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


class OMMGeneral(MinimizationEquilibration):

    @abstractmethod
    def _make_openmm_context(self, system: mm.System) -> mm.Context:
        pass

    @abstractmethod
    def _make_trajectory(self, context: mm.Context) -> mdtraj.Trajectory:
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
            output_context = self._make_openmm_context(system)
            output_traj = self._make_trajectory(output_context)

        return {"trajectory": output_traj, "topology": self.topology}


class OMMMinEq(OMMGeneral):

    def _extract_trajectory(self, context: mm.Context) -> mdtraj.Trajectory:
        final_state = context.getState(getPositions=True, enforcePeriodicBox=False)
        final_coords = final_state.getPositions(asNumpy=True)
        final_vector = final_state.getPeriodicBoxVectors(asNumpy=True)
        output_traj = mdtraj.Trajectory(final_coords[np.newaxis, :] / units.nanometers, topology=self.trajectory.top)
        output_traj.unitcell_vectors = final_vector[np.newaxis, :] / units.nanometers
        return output_traj


class OpenMMMin(OMMMinEq):

    @staticmethod
    def _build_integrator() -> mm.Integrator:
        # return mm.VerletIntegrator(0.1 * units.femtoseconds)
        return mm.LangevinIntegrator(298*units.kelvin,
                                     5/units.picoseconds,
                                     1 * units.femtoseconds)

    def _make_trajectory(self, context: mm.Context) -> mdtraj.Trajectory:
        mm.LocalEnergyMinimizer.minimize(context, tolerance=2)  # , maxIterations=500)
        return self._extract_trajectory(context)

    def _make_openmm_context(self, system: mm.System) -> mm.Context:
        context = mm.Context(system, self._build_integrator(), mm.Platform.getPlatformByName("CPU"))
        context.setPositions(self.trajectory.xyz[-1])
        context.setPeriodicBoxVectors(*self.trajectory.unitcell_vectors[-1])
        return context


class OpenMMEq(OMMMinEq):

    def _make_openmm_context(self, system: mm.System) -> mm.Context:
        integrator = mm.LangevinIntegrator(298*units.kelvin,
                                           5/units.picoseconds,
                                           1 * units.femtoseconds)
        # integrator = mm.VerletIntegrator(1 * units.femtoseconds)
        context = mm.Context(system, integrator)
        context.setPositions(self.trajectory.xyz[-1])
        context.setPeriodicBoxVectors(*self.trajectory.unitcell_vectors[-1])
        context.setVelocitiesToTemperature(298*units.kelvin)
        return context

    def _make_trajectory(self, context: mm.Context) -> mdtraj.Trajectory:
        context.getIntegrator().step(50)
        return self._extract_trajectory(context)

