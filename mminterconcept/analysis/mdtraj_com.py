'''
This module provides access to the MDTraj center of mass function using a MDTraj Trajectory.
'''

from .models import Component
import mdtraj
import numpy

# def mdtraj_com(trajectory_file: str) -> numpy.ndarray:
    # '''
        # This method provides easy access to MDTraj's compute_center_of_mass function.
    # '''
    # trajectory = mdtraj.load_pdb(trajectory_file)
    # return mdtraj.compute_center_of_mass(trajectory)
    
class COMComponent(Component):
    '''
        A component to calculate the Density.
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
        return self.traj.time, mdtraj.compute_center_of_mass(self.traj)
        
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel : str = 'protein'):
        self.process_input(struct, trajectory, sel)
        return self.compute()
