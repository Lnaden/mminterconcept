'''
This module provides access to the a mass density in MDAnalysis using a MDTraj Trajectory.
'''

from .mdanalysis_trajectory_analyzer import MDAnalysisTrajectoryComponent
import numpy
import mdtraj
import MDAnalysis
from .conversion import Distance

from tempfile import TemporaryDirectory

class DensityMDAnalysisComponent(MDAnalysisTrajectoryComponent):
    '''
        A component to calculate the density using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'all'):
        super().__init__(struct, trajectory, str)
    
    def process_input(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel) -> MDAnalysis.Universe:
        return super().process_input(struct, trajectory, sel)
        
    def compute(self) -> numpy.ndarray:
        density_by_frame = numpy.empty(len(self.universe.trajectory))
        time = numpy.ndarray(shape=(len(self.universe.trajectory),))
        mass = self.universe.atoms.select_atoms(self.sel).total_mass()

        for ts in self.universe.trajectory:
            density_by_frame[ts.frame] = mass / ts.volume
            time[ts.frame] = ts.time

        return time, density_by_frame * 1.6605387823355087 / (Distance.ang_to_nm**3)
        
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str='all' ):
        self.process_input(struct, trajectory, sel)
        return self.compute()
