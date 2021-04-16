""" This holds all abstract classes to be implemented elsewhere, for now SingleValueField and Container """

from abc import ABC, abstractmethod

class Container(ABC):
    """ Abstract class to hold atoms in PLUMED. Any Container should:
            -   initialise with a label
            -   return a one-line string for the PLUMED input script
    """
    def __init__(self, label: str):
        self.label = label
        super().__init__()

    @abstractmethod
    def get_plumed_str(self):
        pass


class SingleValueField(ABC):
    """ Abstract class to hold commands in PLUMED which produce a single output. Each SingleValueField should:
            -   initialise with a label
            -   return a one-line string for the PLUMED input script
            -   calculate a value of time-points along the current trajectory
    """
    def __init__(self, label: str):
        self.label = label
        super().__init__()
    
    @abstractmethod
    def get_plumed_str(self):
        pass

    @abstractmethod
    def get_value(self):
        pass