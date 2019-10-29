'''
This module provides access to the MDTraj density function using a MDTraj Trajectory.
'''

from .models import Component
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
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        #self.trajectory = self.process_input(trajectory_file)
        self.trajectory = trajectory
    
    def process_input(self, trajectory: mdtraj.Trajectory) -> mdtraj.Trajectory:
        #return mdtraj.load_pdb(trajectory_file)
        return trajectory
        
    def compute(self):
        return mdtraj.density(self.trajectory)
        # return sum(atom.element.mass for atom in self.trajectory.top.atoms)
        # return self.trajectory.unitcell_volumes
        
    def run(self, trajectory: mdtraj.Trajectory):
        self.process_input(trajectory)
        return self.compute()
        
        
