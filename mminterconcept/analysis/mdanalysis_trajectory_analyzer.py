from mminterconcept.analysis import mda_rmsd, mda_rog, mda_com, mda_density, mda_rdf
import MDAnalysis
import mdtraj

from .models import TrajectoryAnalyzerComponent

from contextlib import contextmanager
import os
from tempfile import TemporaryDirectory

class MDAnalysisTrajectoryComponent(TrajectoryAnalyzerComponent):
    trajectory: mdtraj.Trajectory
    
    def __init__(self, trajectory: mdtraj.Trajectory=None):
        self.trajectory = trajectory
        
    def process_input(self, trajectory: mdtraj.Trajectory) -> mdtraj.Trajectory:
        with TemporaryDirectory() as tempdirname:
            trajectory.save_pdb(tempdirname+'temp.pdb')
            self.universe = MDAnalysis.Universe(tempdirname+'temp.pdb')
        return self.universe
    
    def compute(self):
        pass
    
    def run(self, trajectory: mdtraj.Trajectory, cmd: str):
        cmd_dict = {'rdf': mda_rdf.RDFMDAnalysisComponent, 'rmsd': mda_rmsd.RMSDMDAnalysisComponent, 'rog': mda_rog.ROGMDAnalysisComponent, 'com': mda_com.COMMDAnalysisComponent, 'density': mda_density.DensityMDAnalysisComponent}
        if cmd not in cmd_dict:
            raise KeyError('Operation %s is unsupported by this component.' % (cmd))
        
        comp = cmd_dict[cmd](trajectory)
        return comp.run(trajectory)