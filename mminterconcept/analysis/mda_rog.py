'''
This module provides access to the MDAnalysis radius of gyration function using a MDTraj Trajectory.
'''

from .models import Component
import numpy
import mdtraj
import MDAnalysis
from .conversion import Distance

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory

class ROGMDAnalysisComponent(Component):
    '''
        A component to calculate the radius of gyration using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str='protein'):
        self.trajectory = trajectory
        self.top = top
        self.sel = sel
    
    def process_input(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein') -> MDAnalysis.Universe:
        with TemporaryDirectory() as tempdirname:
            if top:
                top.save(tempdirname+'temp.pdb')
            else:
                trajectory.save(tempdirname+'temp.pdb')

            trajectory.save(tempdirname+'temp.trr')
            self.universe = MDAnalysis.Universe(tempdirname+'temp.pdb', tempdirname+'temp.trr')
        self.sel = sel

        return self.universe
        
    def compute(self):
        rg_by_frame = numpy.empty(len(self.universe.trajectory))
        time = numpy.ndarray(shape=(len(self.universe.trajectory),))

        for ts in self.universe.trajectory:
            rg_by_frame[ts.frame] = self.universe.atoms.select_atoms(self.sel).radius_of_gyration()
            time[ts.frame] = ts.time
        
        return time, rg_by_frame * Distance.ang_to_nm
        
    def run(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein'):
        self.process_input(trajectory, top, sel)
        return self.compute()
