'''
This module provides access to the MDAnalysis rmsd function using a MDTraj Trajectory.
'''

from .models import Component
import numpy as np
import mdtraj
import MDAnalysis
import MDAnalysis.analysis.rms

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory
from .conversion import Distance, Time

class RMSDMDAnalysisComponent(Component):
    '''
        A component to calculate the RMSD using MDAnalysis.
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
            self.universe = MDAnalysis.Universe(tempdirname+'temp.pdb', tempdirname+'temp.trr', superposition=True)
        self.sel = sel

        return self.universe
        
    def compute(self) -> np.ndarray:
        RMSD = MDAnalysis.analysis.rms.RMSD(self.universe.atoms, self.universe, sel=self.sel, groupselections=[self.sel])
        RMSD.run()
        return RMSD.rmsd[:,0] * Time.fs_to_ps, RMSD.rmsd[:,-1] * Distance.ang_to_nm
        
    def run(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str='protein'):
        self.process_input(trajectory, top, sel)
        return self.compute()
