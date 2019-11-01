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
    
    def __init__(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel : str = 'protein'):
        self.trajectory = trajectory
        self.top = top
        self.sel = sel
        
    def process_input(self, trajectory, top: mdtraj.Trajectory = None, sel: str='all'):
        self.traj = trajectory
 
        indices = self.traj.top.select(sel)
        self.traj = self.traj.atom_slice(indices)

        if top:
            self.ref = top.atom_slice(indices)

        return self.traj
        
    def compute(self):
        return self.traj.time, mdtraj.compute_rg(self.traj)
        
    def run(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein'):
        self.process_input(trajectory, top, sel)
        return self.compute()
