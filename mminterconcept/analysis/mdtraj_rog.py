'''
This module provides access to the MDTraj radius of gyration function using a trajectory filepath.
'''

from models import Component
import mdtraj
import numpy

# def mdtraj_rog(trajectory_file: str) -> numpy.ndarray:
    # '''
        # This method provides easy access to MDTraj's compute_rg function.
    # '''
    # trajectory = mdtraj.load_pdb(trajectory_file)
    # return mdtraj.compute_rg(trajectory)
    
class ROGComponent(Component):
    '''
        A component to calculate the radius of gyration.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, trajectory_file: str):
        self.trajectory = self.process_input(trajectory_file)
    
    def process_input(self, trajectory_file: str) -> mdtraj.Trajectory:
        return mdtraj.load_pdb(trajectory_file)
        
    def compute(self):
        return mdtraj.compute_rg(self.trajectory)
        
    def run(self):
        self.process_input(trajectory_file)
        return self.compute()