'''
This module provides access to the MDTraj radius of gyration function using a MDTraj Trajectory.
'''

from .mdtraj_trajectory_analyzer import MDTrajTrajectoryComponent
import mdtraj
import numpy

class ROGComponent(MDTrajTrajectoryComponent):
    '''
        A component to calculate the radius of gyration.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel : str = 'protein'):
        super().__init__(struct, trajectory, sel)

    def process_input(self) -> mdtraj.Trajectory:
        return super().process_input()
        
    def compute(self):
        return self.traj.time, mdtraj.compute_rg(self.traj)
        
    def run(self):
        self.process_input()
        return self.compute()
