import mdtraj

from .models import TrajectoryAnalyzerComponent

class MDTrajTrajectoryComponent(TrajectoryAnalyzerComponent):
    trajectory: mdtraj.Trajectory
    struct: mdtraj.Trajectory

    def __init__(self, struct: mdtraj.Trajectory=None, trajectory: mdtraj.Trajectory=None):
        self.trajectory = trajectory
        self.struct = struct
        
    def process_input(self):
        self.universe = mdtraj.load(self.trajectory, top=self.struct)
        return self.universe

    def compute(self):
        pass
    
    def run(self, trajectory: mdtraj.Trajectory, cmd: str):
        from . import mdtraj_com, mdtraj_rog, mdtraj_rmsd, mdtraj_rdf, mdtraj_density

        cmd_dict = {'rdf': mdtraj_rdf.RDFComponent, 'rmsd': mdtraj_rmsd.RMSDComponent, 'rog': mdtraj_rog.ROGComponent, 'com': mdtraj_com.COMComponent, 'density': mdtraj_density.DensityComponent}
        if cmd not in cmd_dict:
            raise KeyError('Operation %s is unsupported by this component.' % (cmd))
        
        comp = cmd_dict[cmd](trajectory)
        return comp.run(trajectory)
