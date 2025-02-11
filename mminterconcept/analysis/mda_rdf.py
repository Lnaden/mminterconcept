'''
This module provides access to the MDAnalysis RDF function using a MDTraj Trajectory.
'''
from .models import Component
import numpy
import mdtraj
import MDAnalysis
from MDAnalysis.analysis.rdf import InterRDF

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory
from .conversion import Distance

class RDFMDAnalysisComponent(Component):
    '''
        A component to calculate the density using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein'):
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
        group = self.universe.select_atoms(self.sel)
        rdf = InterRDF(group, group, nbins=200, range=[0.0, 10.0])
        rdf.run()
        bins = rdf.bins
        return bins * Distance.ang_to_nm, rdf.rdf

    def run(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'protein'):
        self.process_input(trajectory, top, sel)
        return self.compute()
