'''
This module provides access to the MDAnalysis radius of gyration function using a MDTraj Trajectory.
'''

from .mdanalysis_trajectory_analyzer import MDAnalysisTrajectoryComponent
import numpy
import mdtraj
import MDAnalysis
from .conversion import Distance

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory

class ROGMDAnalysisComponent(MDAnalysisTrajectoryComponent):
    '''
        A component to calculate the radius of gyration using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str='protein'):
        super().__init__(struct, trajectory, sel)
    
    def process_input(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str='protein') -> mdtraj.Trajectory:
        return super().process_input(struct, trajectory, sel)
        
    def compute(self):
        rg_by_frame = numpy.empty(len(self.universe.trajectory))
        time = numpy.ndarray(shape=(len(self.universe.trajectory),))

        for ts in self.universe.trajectory:
            rg_by_frame[ts.frame] = self.universe.atoms.select_atoms(self.sel).radius_of_gyration()
            time[ts.frame] = ts.time
        
        return time, rg_by_frame * Distance.ang_to_nm
        
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        self.process_input(struct, trajectory, sel)
        return self.compute()
