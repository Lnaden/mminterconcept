'''
This module provides access to the MDTraj center of mass function using a MDTraj Trajectory.
'''

from .mdtraj_trajectory_analyzer import MDTrajTrajectoryComponent
import mdtraj
import numpy

# def mdtraj_com(trajectory_file: str) -> numpy.ndarray:
    # '''
        # This method provides easy access to MDTraj's compute_center_of_mass function.
    # '''
    # trajectory = mdtraj.load_pdb(trajectory_file)
    # return mdtraj.compute_center_of_mass(trajectory)
    
class COMComponent(MDTrajTrajectoryComponent):
    '''
        A component to calculate the Density.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory):
        super().__init__(struct, trajectory)
 
    def process_input(self) -> mdtraj.Trajectory:
        return super().process_input(self.struct, self.trajectory)
        
    def compute(self):
        return mdtraj.compute_center_of_mass(self.universe)
        
    def run(self):
        self.process_input()
        return self.compute()
