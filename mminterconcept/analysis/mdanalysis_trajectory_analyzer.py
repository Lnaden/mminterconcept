import MDAnalysis
import mdtraj

from .models import TrajectoryAnalyzerComponent
from tempfile import TemporaryDirectory


class MDAnalysisTrajectoryComponent(TrajectoryAnalyzerComponent):
    trajectory: mdtraj.Trajectory
    
    def __init__(self, trajectory: mdtraj.Trajectory):
        self.trajectory = trajectory

    def process_input(self) -> MDAnalysis.Universe:
        with TemporaryDirectory() as tempdirname:
            self.trajectory.save(tempdirname+'temp.pdb')
            self.trajectory.save(tempdirname+'temp.trr')
            self.universe = MDAnalysis.Universe(tempdirname+'temp.pdb', tempdirname+'temp.trr')

        return self.universe
    
    def compute(self):
        raise NotImplementedError
    
    def run(self, cmd: str):
        from . import mda_rmsd, mda_rog, mda_com, mda_density, mda_rdf

        cmd_dict = {'rdf': mda_rdf.RDFMDAnalysisComponent, 'rmsd': mda_rmsd.RMSDMDAnalysisComponent, 'rog': mda_rog.ROGMDAnalysisComponent, 
                    'com': mda_com.COMMDAnalysisComponent, 'density': mda_density.DensityMDAnalysisComponent}

        if cmd not in cmd_dict:
            raise KeyError(f'Operation {cmd} not supported by this component.')
        
        comp = cmd_dict[cmd]
        return comp.run()
