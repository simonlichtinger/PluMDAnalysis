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



class PLUMED_Combine(SingleValueField):
    def __init__(self, arguments:list, coefficients:list, label:str):
        """Constructure for the PLUMED_Combine class. Initialised with argument SingleValueFields and the coefficients.

        :param arguments:       List of SingleValueFields to be used as arguments for the linear combination.
        :type arguments:        List<SingleValueField>
        :param coefficients:    List of coefficients for each 
        :type coefficients:     List<float>
        :param label:           Label for this PLUMED_Combine.
        :type label:            str
        """
        super(PLUMED_Combine, self).__init__(label)
        if not len(arguments) == len(coefficients):
            raise ValueError("Arguments and Coefficients need to have the same length.")
        self.arguments, self.coefficients = arguments, coefficients

    def get_plumed_str(self):
        """Obtain the PLUMED representation of this PLUMED_Combine as a string.

        :return:    The string representing this PLUMED_Combine, ready to be placed in input script.
        :rtype:     str
        """
        return f"{self.label}: COMBINE ARG={','.join([x.label for x in self.arguments])} COEFFICIENTS={','.join([str(x) for x in self.coefficients])} PERIODIC=NO"

    def get_value(self):
        """Obtain the instantaneous value of this linear combination for the current time frame.

        :return:    Value of the linear combination.
        :rtype:     float
        """
        return np.sum(np.array([arg.get_value() for arg in self.arguments]).dot(self.coefficients))
