'''
This module provides access to the a mass density in MDAnalysis using a MDTraj Trajectory.
'''

from .models import Component
import numpy
import mdtraj
import MDAnalysis
from .conversion import Distance

from tempfile import TemporaryDirectory

class DensityMDAnalysisComponent(Component):
    '''
        A component to calculate the density using MDAnalysis.
    '''
    
    universe: MDAnalysis.Universe
    
    def __init__(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'all'):
        self.trajectory = trajectory
        self.top = top
        self.sel = sel
    
    def process_input(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str = 'all') -> MDAnalysis.Universe:
        with TemporaryDirectory() as tempdirname:

            if top:
                top.save(tempdirname+'temp.pdb')
            else:
                trajectory.save(tempdirname+'temp.pdb')

            trajectory.save(tempdirname+'temp.trr')
            self.universe = MDAnalysis.Universe(tempdirname+'temp.pdb', tempdirname+'temp.trr')
        self.sel = sel

        return self.universe
        
    def compute(self) -> numpy.ndarray:
        density_by_frame = numpy.empty(len(self.universe.trajectory))
        time = numpy.ndarray(shape=(len(self.universe.trajectory),))
        mass = self.universe.atoms.select_atoms(self.sel).total_mass()

        for ts in self.universe.trajectory:
            density_by_frame[ts.frame] = mass / ts.volume
            time[ts.frame] = ts.time

        return time, density_by_frame * 1.6605387823355087 / (Distance.ang_to_nm**3)
        
    def run(self, trajectory: mdtraj.Trajectory, top: mdtraj.Trajectory = None, sel: str='all' ):
        self.process_input(trajectory, top, sel)
        return self.compute()
