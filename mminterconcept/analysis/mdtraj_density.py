'''
This module provides access to the MDTraj density function using a MDTraj Trajectory.
'''

from .mdtraj_trajectory_analyzer import MDTrajTrajectoryComponent
import mdtraj
import numpy
from typing import Tuple

# def mdtraj_density(trajectory_file: str) -> numpy.ndarray:
    # '''
        # This method provides easy access to MDTraj's density function.
    # '''
    # trajectory = mdtraj.load_pdb(trajectory_file)
    # return mdtraj.density(trajectory)
    
class DensityComponent(MDTrajTrajectoryComponent):
    '''
        A component to calculate the Density.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'all'):
        super().__init__(struct, trajectory, sel)
    
    def process_input(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str) -> mdtraj.Trajectory:
        return super().process_input(struct, trajectory, sel)
        
    def compute(self) -> Tuple[numpy.ndarray, numpy.ndarray]:
        masses = numpy.array([atom.element.mass for atom in self.traj.top.atoms])
        return self.traj.time, mdtraj.density(self.traj, masses)
        
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str='all'):
        self.process_input(struct, trajectory, sel)
        return self.compute()
