'''
This module provides access to the MDAnalysis radius of gyration function using a MDTraj Trajectory.
'''

from models import Component
import numpy
import mdtraj
import MDAnalysis

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory

class ROGMDAnalysisComponent(Component):
    '''
        A component to calculate the radius of gyration using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        self.universe = self.process_input(trajectory)
    
    def process_input(self, trajectory: mdtraj.Trajectory) -> mdtraj.Trajectory:
        with TemporaryDirectory() as tempdirname:
            trajectory.save_pdb(tempdirname+'temp.pdb')
            self.universe = MDAnalysis.Universe(tempdirname+'temp.pdb')
        return self.universe
        
    def compute(self):
        rg_by_frame = numpy.empty(len(self.universe.trajectory))
        for ts in self.universe.trajectory:
            rg_by_frame[ts.frame] = self.universe.atoms.radius_of_gyration() / 10.0
        return rg_by_frame
        
    def run(self, trajectory: mdtraj.Trajectory):
        self.process_input(trajectory)
        return self.compute()