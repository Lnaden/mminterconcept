'''
This module provides access to the MDTraj radius of gyration function using a MDTraj Trajectory.
'''

from .models import Component
import mdtraj
import numpy

class ROGComponent(Component):
    '''
        A component to calculate the radius of gyration.
    '''
    trajectory: mdtraj.Trajectory
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel : str = 'protein'):
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
        return self.traj.time, mdtraj.compute_rg(self.traj)
        
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        self.process_input(struct, trajectory, sel)
        return self.compute()
