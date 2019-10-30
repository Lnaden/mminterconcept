'''
This module provides access to the MDAnalysis RDF function using a MDTraj Trajectory.
'''
from .mdanalysis_trajectory_analyzer import MDAnalysisTrajectoryComponent
import numpy
import mdtraj
import MDAnalysis
from MDAnalysis.analysis.rdf import InterRDF

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory
from .conversion import Distance

class RDFMDAnalysisComponent(MDAnalysisTrajectoryComponent):
    '''
        A component to calculate the density using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        super().__init__(struct, trajectory, sel)
    
    def process_input(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str) -> mdtraj.Trajectory:
        return super().process_input(struct, trajectory, sel)
        
    def compute(self):
        group = self.universe.select_atoms(self.sel)
        rdf = InterRDF(group, group, nbins=200, range=[0.0, 10.0])
        rdf.run()
        bins = rdf.bins
        return bins * Distance.ang_to_nm, rdf.rdf

    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str = 'protein'):
        self.process_input(struct, trajectory, sel)
        return self.compute()
