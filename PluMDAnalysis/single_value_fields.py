""" This holds all implementations of the SingleValueField parent class """

import numpy as np
from .containers import PLUMED_Group
from .specifications import SingleValueField

class PLUMED_Distance(SingleValueField):
    """ Class to hold two AtomGroup objects, calculate geometric distances between them
     and transform the distance into a PLUMED input string, implements SingleValueField. """
    def __init__(self, group_1: PLUMED_Group, group_2: PLUMED_Group, label: str):
        """Constructor for the PLUMED_Distance class, initialised with two PLUMED_Group objects.

        :param group_1:         First group of two defining the distance.
        :type group_1:          PLUMED_Group
        :param group_2:         Second group of two defining the distance.
        :type group_2:          PLUMED_Group
        :param label:           Label for this PLUMED_Distance.
        :type label:            str
        """
        self.group_1, self.group_2 = group_1, group_2
        super(PLUMED_Distance, self).__init__(label)

    def get_plumed_str(self):
        """Obtain the PLUMED representation of this PLUMED_Distance as a string.

        :return: The string representing this PLUMED_Distance, ready to be placed in input script.
        :rtype: str
        """
        return self.label + ": DISTANCE ATOMS=" + self.group_1.label + "," + self.group_2.label

    def get_value(self):
        """Obtain the geometric distance between the two groups, at the current frame, in nm.

        :return:    Distance between the two groups, at the current frame, in nm.
        :rtype:     float
        """
        return np.sqrt(np.sum((self.group_1.atom_group.center_of_geometry() - self.group_2.atom_group.center_of_geometry())**2)) / 10