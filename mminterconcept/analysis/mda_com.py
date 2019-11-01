'''
This module provides access to the MDAnalysis center of mass function using a MDTraj Trajectory.
'''

from .models import Component
import numpy
import mdtraj
import MDAnalysis

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory

from .conversion import Distance, Time

class COMMDAnalysisComponent(Component):
    '''
        A component to calculate the center of mass using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein'):
        self.trajectory = trajectory
        self.top = top
        self.sel = sel

    def process_input(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein') -> MDAnalysis.Universe:
        with TemporaryDirectory() as tempdirname:

            if top:
                struct.save(tempdirname+'temp.pdb')
            else:
                trajectory.save(tempdirname+'temp.pdb')

            trajectory.save(tempdirname+'temp.trr')

            self.universe = MDAnalysis.Universe(tempdirname+'temp.pdb', tempdirname+'temp.trr')
        self.sel = sel

        return self.universe

    def compute(self) -> numpy.ndarray:
        com_by_frame = numpy.ndarray(shape=(len(self.universe.trajectory), numpy.size(self.universe.atoms.positions,1)))
        time = numpy.ndarray(shape=(len(self.universe.trajectory),))

        group = self.universe.atoms.select_atoms(self.sel)

        for ts in self.universe.trajectory:
            com_by_frame[ts.frame,:] = group.center_of_mass()
            time[ts.frame] = ts.time

        return time, com_by_frame * Distance.ang_to_nm
        
    def run(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory, sel: str='protein') -> numpy.ndarray:
        self.process_input(trajectory, top, sel)
        return self.compute()
