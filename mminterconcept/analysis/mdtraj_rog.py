'''
This module provides access to the MDTraj radius of gyration function using a MDTraj Trajectory.
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
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        #self.trajectory = self.process_input(trajectory_file)
        self.trajectory = trajectory
    
    def process_input(self, trajectory: mdtraj.Trajectory) -> mdtraj.Trajectory:
        #return mdtraj.load_pdb(trajectory_file)
        return trajectory
        
    def compute(self):
        return mdtraj.compute_rg(self.trajectory)
        
    def run(self, trajectory: mdtraj.Trajectory):
        self.process_input(trajectory)
        return self.compute()