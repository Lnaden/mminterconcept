'''
This module provides access to the MDTraj density function using a MDTraj Trajectory.
'''

from .models import Component
import mdtraj
import numpy
from typing import Tuple

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
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'all'):
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
        
    def compute(self) -> Tuple[numpy.ndarray, numpy.ndarray]:
        masses = numpy.array([atom.element.mass for atom in self.traj.top.atoms])
        return self.traj.time, mdtraj.density(self.traj, masses)
        
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str='all'):
        self.process_input(struct, trajectory, sel)
        return self.compute()
