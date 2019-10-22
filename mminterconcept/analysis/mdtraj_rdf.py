'''
This module provides access to the MDTraj rdf function using a trajectory filepath.
'''

from models import Component
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
    
    def __init__(self, trajectory_file: str):
        self.trajectory = self.process_input(trajectory_file)
    
    def process_input(self, trajectory_file: str) -> mdtraj.Trajectory:
        return mdtraj.load_pdb(trajectory_file)
        
    def compute(self):
        pairs = self.trajectory.top.select_pairs('all', 'all')
        return mdtraj.compute_rdf(self.trajectory, pairs)
        
    def run(self):
        self.process_input(trajectory_file)
        return self.compute()