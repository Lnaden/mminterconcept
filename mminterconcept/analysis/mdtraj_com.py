'''
This module provides access to the MDTraj center of mass function using a trajectory filepath.
'''

from models import Component
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
    
    def __init__(self, trajectory_file: str):
        self.trajectory = self.process_input(trajectory_file)
    
    def process_input(self, trajectory_file: str) -> mdtraj.Trajectory:
        return mdtraj.load_pdb(trajectory_file)
        
    def compute(self):
        return mdtraj.compute_center_of_mass(self.trajectory)
        
    def run(self):
        self.process_input(trajectory_file)
        return self.compute()