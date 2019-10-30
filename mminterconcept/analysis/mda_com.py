'''
This module provides access to the MDAnalysis center of mass function using a MDTraj Trajectory.
'''

from .mdanalysis_trajectory_analyzer import MDAnalysisTrajectoryComponent
import numpy
import mdtraj
import MDAnalysis

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory

from .conversion import Distance, Time

class COMMDAnalysisComponent(MDAnalysisTrajectoryComponent):
    '''
        A component to calculate the center of mass using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        super().__init__(struct, trajectory, sel)
    
    def process_input(self, struct, trajectory, sel) -> MDAnalysis.Universe:
        return super().process_input(struct, trajectory, sel)

    def compute(self) -> numpy.ndarray:
        com_by_frame = numpy.ndarray(shape=(len(self.universe.trajectory), numpy.size(self.universe.atoms.positions,1)))
        time = numpy.ndarray(shape=(len(self.universe.trajectory),))

        group = self.universe.atoms.select_atoms(self.sel)

        for ts in self.universe.trajectory:
            com_by_frame[ts.frame,:] = group.center_of_mass()
            time[ts.frame] = ts.time

        return time, com_by_frame * Distance.ang_to_nm
        
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str='protein') -> numpy.ndarray:
        self.process_input(struct, trajectory, sel)
        return self.compute()
