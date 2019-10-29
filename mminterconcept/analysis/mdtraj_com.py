'''
This module provides access to the MDTraj center of mass function using a MDTraj Trajectory.
'''

from .models import Component
import mdtraj
import numpy

# def mdtraj_com(trajectory_file: str) -> numpy.ndarray:
    # '''
        # This method provides easy access to MDTraj's compute_center_of_mass function.
    # '''
    # trajectory = mdtraj.load_pdb(trajectory_file)
    # return mdtraj.compute_center_of_mass(trajectory)
    
class COMComponent(Component):
    '''
        A component to calculate the Density.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        #self.trajectory = self.process_input(trajectory_file)
        self.trajectory = trajectory
    
    def process_input(self, trajectory: mdtraj.Trajectory) -> mdtraj.Trajectory:
        #return mdtraj.load_pdb(trajectory_file)
        return trajectory
        
    def compute(self):
        return mdtraj.compute_center_of_mass(self.trajectory)
        
    def run(self, trajectory: mdtraj.Trajectory):
        self.process_input(trajectory)
        return self.compute()
