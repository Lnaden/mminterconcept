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
    
    def __init__(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein'):
        self.trajectory = trajectory
        self.top = top
        self.sel = sel
        
    def process_input(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel='all'):
        self.traj = trajectory
 
        indices = self.traj.top.select(sel)
        self.traj = self.traj.atom_slice(indices)

        if top:
            self.ref = top.atom_slice(indices)

        return self.traj

    def compute(self):
        pairs = self.traj.top.select_pairs('all', 'all')
        return mdtraj.compute_rdf(self.traj, pairs)

    def run(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein'):
        self.process_input(trajectory, top, sel)
        return self.compute()
