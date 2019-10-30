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

class COMMDAnalysisComponent(MDAnalysisTrajectoryComponent):
    '''
        A component to calculate the center of mass using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        super().__init__(trajectory)
    
    def process_input(self) -> MDAnalysis.Universe:
        return super().process_input()

    def compute(self) -> numpy.ndarray:
        com_by_frame = numpy.ndarray(shape=(len(self.universe.trajectory), numpy.size(self.universe.atoms.positions,1)))
        for ts in self.universe.trajectory:
            com_by_frame[ts.frame,:] = self.universe.atoms.center_of_mass() / 10.0
        return com_by_frame
        
    def run(self) -> numpy.ndarray:
        self.process_input()
        return self.compute()
