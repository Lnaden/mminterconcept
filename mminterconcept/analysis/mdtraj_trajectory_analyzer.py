from mminterconcept.analysis import mdtraj_com, mdtraj_rog, mdtraj_rmsd, mdtraj_rdf, mdtraj_density
import mdtraj

from .models import TrajectoryAnalyzerComponent

class MDTrajTrajectoryComponent(TrajectoryAnalyzerComponent):
    trajectory: mdtraj.Trajectory
    
    def __init__(self, trajectory: mdtraj.Trajectory=None):
        self.trajectory = trajectory
        
    def process_input(self, trajectory: mdtraj.Trajectory):
        self.trajectory = trajectory
        return self.trajectory
    
    def compute(self):
        pass
    
    def run(self, trajectory: mdtraj.Trajectory, cmd: str):
        cmd_dict = {'rdf': mdtraj_rdf.RDFComponent, 'rmsd': mdtraj_rmsd.RMSDComponent, 'rog': mdtraj_rog.ROGComponent, 'com': mdtraj_com.COMComponent, 'density': mdtraj_density.DensityComponent}
        if cmd not in cmd_dict:
            raise KeyError('Operation %s is unsupported by this component.' % (cmd))
        
        comp = cmd_dict[cmd](trajectory)
        return comp.run(trajectory)