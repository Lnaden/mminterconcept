'''
This module provides access to the MDAnalysis rmsd function using a MDTraj Trajectory.
'''

from .models import Component
import numpy
import mdtraj
import MDAnalysis
import MDAnalysis.analysis.rms

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory

# @contextmanager
# def tempcd(*args, **kwargs):
    # cur = os.getcwd()
    # try:
        # with TemporaryDirectory(*args, **kwargs) as tmpdir:
            # os.chdir(tmpdir)
            # yield
    # finally:
        # os.chdir(cur)
            
class RMSDMDAnalysisComponent(Component):
    '''
        A component to calculate the RMSD using MDAnalysis.
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
        protein = self.universe.select_atoms('protein')
        a = protein.positions.copy()
        self.universe.trajectory[-1]
        b = protein.positions.copy()
        rmsd = MDAnalysis.analysis.rms.rmsd(a, b, center=True)
        return rmsd
        
    def run(self, trajectory: mdtraj.Trajectory):
        self.process_input(trajectory)
        return self.compute()
