'''
This module provides access to the a mass density in MDAnalysis using a MDTraj Trajectory.
'''

from .models import Component
import numpy
import mdtraj
import MDAnalysis

from tempfile import TemporaryDirectory

class DensityMDAnalysisComponent(Component):
    '''
        A component to calculate the density using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        super().__init__(trajectory)
    
    def process_input(self, trajectory: mdtraj.Trajectory) -> MDAnalysis.Universe:
        return super().process_input(trajectory)
        
    def compute(self) -> numpy.ndarray:
        density_by_frame = numpy.empty(len(self.universe.trajectory))
        mass = self.universe.atoms.total_mass()
        for ts in self.universe.trajectory:
            density_by_frame[ts.frame] = (mass / (ts.volume / 1000.)) * 1.6605387823355087
        return density_by_frame
        
    def run(self, trajectory: mdtraj.Trajectory):
        self.process_input(trajectory)
        return self.compute()
