'''
This module provides a base component API to be extended by components.
'''

from abc import ABC, abstractmethod

class Component(ABC):
    def __init__():
        pass
        
    @abstractmethod
    def process_input():
        pass
    
    @abstractmethod
    def compute():
        pass
    
    @abstractmethod
    def run():
        pass