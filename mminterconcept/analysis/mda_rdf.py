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

class RDFMDAnalysisComponent(Component):
    '''
        A component to calculate the density using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        super().__init__(trajectory)
    
    def process_input(self, trajectory: mdtraj.Trajectory) -> mdtraj.Trajectory:
        return super().process_input(trajectory)
        
    def compute(self):
        # density_by_frame = numpy.empty(len(self.universe.trajectory))
        # mass = self.universe.atoms.total_mass()
        # for ts in self.universe.trajectory:
            # density_by_frame[ts.frame] = (mass / (ts.volume / 1000)) * 1.6605387823355087
        # return density_by_frame
        group = self.universe.select_atoms('protein')
        # ags = [[group, group]]
        # print(len(group))
        rdf = InterRDF(group, group, nbins=200, range=[0.0, 10.0])
        rdf.run()
        bins = rdf.bins
        return bins, rdf.rdf

    def run(self, trajectory: mdtraj.Trajectory):
        self.process_input(trajectory)
        return self.compute()
