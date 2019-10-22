'''
This module provides access to the MDTraj rmsd function using a trajectory filepath.
'''

from models import Component
import numpy
import mdtraj

# def mdtraj_rmsd(trajectory_file: str) -> numpy.ndarray:
    # '''
        # This method provides easy access to MDTraj's rmsd function.
    # '''
    # trajectory = mdtraj.load_pdb(trajectory_file)
    # return mdtraj.rmsd(trajectory, trajectory, 0)
    
class RMSDComponent(Component):
    '''
        A component to calculate the RMSD.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, trajectory_file: str):
        self.trajectory = self.process_input(trajectory_file)
    
    def process_input(self, trajectory_file: str) -> mdtraj.Trajectory:
        return mdtraj.load_pdb(trajectory_file)
        
    def compute(self):
        return mdtraj.rmsd(self.trajectory, self.trajectory, 0)
        
    def run(self):
        self.process_input(trajectory_file)
        return self.compute()