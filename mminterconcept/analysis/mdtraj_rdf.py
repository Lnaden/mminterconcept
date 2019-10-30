'''
This module provides access to the MDTraj rdf function using a MDTraj Trajectory.
'''

from .mdtraj_trajectory_analyzer import MDTrajTrajectoryComponent
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
    
    
class RDFComponent(MDTrajTrajectoryComponent):
    '''
        A component to calculate the RDF.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        super().__init__(struct, trajectory, sel)
    
    def process_input(self) -> mdtraj.Trajectory:
        super().process_input()

    def compute(self):
        pairs = self.traj.top.select_pairs('all', 'all')
        return mdtraj.compute_rdf(self.traj, pairs)

    def run(self):
        self.process_input()
        return self.compute()
