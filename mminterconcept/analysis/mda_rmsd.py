'''
This module provides access to the MDAnalysis rmsd function using a MDTraj Trajectory.
'''

from .mdanalysis_trajectory_analyzer import MDAnalysisTrajectoryComponent
import numpy as np
import mdtraj
import MDAnalysis
import MDAnalysis.analysis.rms

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory

# @contextmanager
# def tempcd(*args, **kwargs):
    # cur = os.getcwd()
    # try:
        # with TemporaryDirectory(*args, **kwargs) as tmpdir:
            # os.chdir(tmpdir)
            # yield
    # finally:
        # os.chdir(cur)
            
class RMSDMDAnalysisComponent(MDAnalysisTrajectoryComponent):
    '''
        A component to calculate the RMSD using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        super().__init__(trajectory)
    
    def process_input(self) -> MDAnalysis.Universe:
        return super().process_input(self.trajectory)
        
    def compute(self) -> np.ndarray:
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe.atoms, self.universe, sel='protein')
        RMSD.run()
        return RMSD.rmsd.flatten()
        
    def run(self):
        self.process_input()
        return self.compute()
