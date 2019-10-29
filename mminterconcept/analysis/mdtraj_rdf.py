'''
This module provides access to the MDTraj rdf function using a MDTraj Trajectory.
'''

from .models import Component
import numpy
import mdtraj
from typing import Tuple

# def mdtraj_rdf(trajectory_file: str) -> Tuple[numpy.ndarray, numpy.ndarray]:
    # '''
        # This method provides easy access to MDTraj's compute_rdf function.
    # '''
    # trajectory = mdtraj.load_pdb(trajectory_file)
    # pairs = trajectory.top.select_pairs('all', 'all')
    # return mdtraj.compute_rdf(trajectory, pairs)
    
    
class RDFComponent(Component):
    '''
        A component to calculate the RDF.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        #self.trajectory = self.process_input(trajectory_file)
        self.trajectory = trajectory
    
    def process_input(self, trajectory: mdtraj.Trajectory) -> mdtraj.Trajectory:
        #return mdtraj.load_pdb(trajectory_file)
        return trajectory
        
    def compute(self):
        #water = self.trajectory.top.select('water')
        pairs = self.trajectory.top.select_pairs('name C', 'name C')
        return mdtraj.compute_rdf(self.trajectory, pairs)
        
    def run(self, trajectory: mdtraj.Trajectory):
        self.process_input(trajectory)
        return self.compute()
