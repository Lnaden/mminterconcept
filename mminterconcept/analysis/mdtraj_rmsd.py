'''
This module provides access to the MDTraj rmsd function using a MDTraj Trajectory.
'''

from .mdtraj_trajectory_analyzer import MDTrajTrajectoryComponent
import numpy
import mdtraj

class RMSDComponent(MDTrajTrajectoryComponent):
    '''
        A component to calculate the RMSD.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        super().__init__(struct, trajectory, sel)
    
    def process_input(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str) -> mdtraj.Trajectory:
        return super().process_input(struct, trajectory, sel)

    def compute(self):
        return self.traj.time, mdtraj.rmsd(self.traj, self.ref)
        
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        self.process_input(struct, trajectory, sel)
        return self.compute()
