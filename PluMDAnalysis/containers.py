""" This holds all implementations of the Container parent class. """

import MDAnalysis
from MDAnalysis.core.groups import AtomGroup
from .specifications import Container

class PLUMED_Group(Container):
    """ Class to hold an AtomGroup and transform it into a PLUMED input string, implements Container. """
    def __init__(self, atom_group: AtomGroup, label: str, group_type='CENTER'):
        """Construcutor for the PLUMED_Group class, initialised with an AtomGroup.

        :param atom_group:  AtomGroup to be held by this instance of PLUMED_Group.
        :type atom_group:   AtomGroup
        :param label:       Label for this instance of PLUMED_Group.
        :type label:        str
        :param group_type:  Type of PLUMED group, 'CENTER' or 'COM', defaults to 'CENTER'
        :type group_type:   str, optional
        :raises TypeError:  The provided atom_group parameter is not of type AtomGroup.
        :raises valueError:  The provided atom_group parameter must be 'CENTER' or 'COM'.
        """
        super(PLUMED_Group, self).__init__(label)
        if not isinstance(atom_group, AtomGroup):
            raise TypeError("The provided atom_group parameter is not of type AtomGroup.")
        self.atom_group = atom_group
        if not group_type in ['CENTER', 'COM']:
            raise ValueError("The provided atom_group parameter must be 'CENTER' or 'COM'.")
        self.group_type = group_type

    def get_plumed_str(self):
        """Obtain the PLUMED representation of this PLUMED_Group as a string.

        :return: The string representing this PLUMED_Group, ready to be placed in input script.
        :rtype: str
        """
        output_str = self.label + ": " + self.group_type + " ATOMS=" + ','.join([str(atom.ix+1) for atom in self.atom_group])
        return output_str
