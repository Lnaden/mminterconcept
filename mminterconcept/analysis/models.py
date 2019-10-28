'''
This module provides a base component API to be extended by components.
'''

from abc import ABC, abstractmethod
import mdtraj

class Component(ABC):
    def __init__():
        raise NotImplementedError('__init__ for Abstract class for a Component.')
        
    @abstractmethod
    def process_input():
        raise NotImplementedError('process_input for Abstract class for a Component.')
    
    @abstractmethod
    def compute():
        raise NotImplementedError('compute for Abstract class for a Component.')
    
    @abstractmethod
    def run():
        raise NotImplementedError('run for Abstract class for a Component.')
        
class TrajectoryAnalyzerComponent(Component):

    @abstractmethod
    def __init__(self, trajectory: mdtraj.Trajectory=None):
        raise NotImplementedError('Abstract base class for Trajectory Analysis.')
    
    @abstractmethod
    def process_input(self, trajectory: mdtraj.Trajectory):
        raise NotImplementedError('Abstract base class for Trajectory Analysis.')
        
    @abstractmethod
    def compute(self):
        raise NotImplementedError('Abstract base class for Trajectory Analysis.')
    
    @abstractmethod
    def run(self, trajectory: mdtraj.Trajectory, cmd: str):
        raise NotImplementedError('Abstract base class for Trajectory Analysis.')