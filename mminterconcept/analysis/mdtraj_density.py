'''
This module provides access to the MDTraj density function using a trajectory filepath.
'''

from models import Component
import mdtraj
import numpy

# def mdtraj_density(trajectory_file: str) -> numpy.ndarray:
    # '''
        # This method provides easy access to MDTraj's density function.
    # '''
    # trajectory = mdtraj.load_pdb(trajectory_file)
    # return mdtraj.density(trajectory)
    
class DensityComponent(Component):
    '''
        A component to calculate the Density.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, trajectory_file: str):
        self.trajectory = self.process_input(trajectory_file)
    
    def process_input(self, trajectory_file: str) -> mdtraj.Trajectory:
        return mdtraj.load_pdb(trajectory_file)
        
    def compute(self):
        return mdtraj.density(self.trajectory)
        
    def run(self):
        self.process_input(trajectory_file)
        return self.compute()
        
        