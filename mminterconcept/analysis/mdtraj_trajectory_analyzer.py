import mdtraj

from .models import TrajectoryAnalyzerComponent

class MDTrajTrajectoryComponent(TrajectoryAnalyzerComponent):
    trajectory: mdtraj.Trajectory
    struct: mdtraj.Trajectory

    def __init__(self, struct: mdtraj.Trajectory=None, trajectory: mdtraj.Trajectory=None, sel : str = 'all'):
        self.trajectory = trajectory
        self.struct = struct
        self.sel = sel
        
    def process_input(self, struct, trajectory, sel='all'):
        self.traj = trajectory
        self.ref = struct
 
        indices = self.traj.top.select(sel)
        self.traj = self.traj.atom_slice(indices)
        self.ref = self.ref.atom_slice(indices)

        return self.traj

    def compute(self):
        pass
    
    def run(self, struct: mdtraj.Trajectory, trajectory: mdtraj.Trajectory, sel: str, cmd: str):
        from . import mdtraj_com, mdtraj_rog, mdtraj_rmsd, mdtraj_rdf, mdtraj_density

        cmd_dict = {'rdf': mdtraj_rdf.RDFComponent, 'rmsd': mdtraj_rmsd.RMSDComponent, 'rog': mdtraj_rog.ROGComponent, 'com': mdtraj_com.COMComponent, 'density': mdtraj_density.DensityComponent}
        if cmd not in cmd_dict:
            raise KeyError(f'Operation {cmd} is unsupported by this component.')
        
        comp = cmd_dict[cmd]
        return comp.run(struct, trajectory, sel)
