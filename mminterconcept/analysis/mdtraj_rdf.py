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
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        self.trajectory = trajectory
        self.struct = struct
        self.sel = sel
        
    def process_input(self, struct, trajectory, sel='all'):
        self.traj = trajectory
        self.ref = struct
 
        indices = self.traj.top.select(sel)
        self.traj = self.traj.atom_slice(indices)
        self.ref = self.ref.atom_slice(indices)

        return self.traj

    def compute(self):
        pairs = self.traj.top.select_pairs('all', 'all')
        return mdtraj.compute_rdf(self.traj, pairs)

    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        self.process_input(struct, trajectory, sel)
        return self.compute()
