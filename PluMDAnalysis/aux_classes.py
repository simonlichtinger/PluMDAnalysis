""" This holds the code for auxiliary classes, including PLUMED_group, Distance and Restraint. """

import MDAnalysis
from MDAnalysis.core.groups import AtomGroup

import numpy as np

class PLUMED_Group:
    """ Class to hold an AtomGroup and transform it into a PLUMED input string. """
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
        if not isinstance(atom_group, AtomGroup):
            raise TypeError("The provided atom_group parameter is not of type AtomGroup.")
        self.atom_group = atom_group
        self.label = label
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


class PLUMED_Distance:
    """ Class to hold two AtomGroup objects, calculate geometric distances between them
     and transform the distance into a PLUMED input string. """
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
        self.label = label

    def get_plumed_str(self):
        """Obtain the PLUMED representation of this PLUMED_Distance as a string.

        :return: The string representing this PLUMED_Distance, ready to be placed in input script.
        :rtype: str
        """
        return self.label + ": DISTANCE ATOMS=" + self.group_1.label + "," + self.group_2.label

    def get_distance(self):
        """Obtain the geometric distance between the two groups, at the current frame, in nm.

        :return:    Distance between the two groups, at the current frame, in nm.
        :rtype:     float
        """
        return np.sqrt(np.sum((self.group_1.atom_group.center_of_geometry() - self.group_2.atom_group.center_of_geometry())**2)) / 10


class PLUMED_Restraint:
    """ Class to define a PLUMED moving restraint, and output it as a list of strings. """
    def __init__(self, distances: list):
        """Constructor with the PLUMED_Distance objects to be included in the restraint.

        :param distances:   List of PLUMED_Distance objects to be included in the restraint.
        :type distances:    List<PLUMED_Distance>
        """
        self.distances = distances

    def add_time_series(self, time_series: list, kappas: list):
        """Adds data for time-steps and kappa values across which the restraint is to work.

        :param time_series:  List of time-series points at which the PLUMED moving restraint takes steps.
        :type time_series:   List<int>
        :param kappas:       List of lists to record all the kappa values required for the restraint. First dimension are time-steps, second dimension are the distances.
        :type kappas:        List<List<int>>
        :raises ValueError:  Inconsistent input for restraint generation (some consistency checks)
        """
        # Do some consistency checks
        if not len(time_series) == len(kappas):
            raise ValueError("Inconsistent input for restraint generation: lenghts of time_series and kappa do not match.")
        for step in kappas:
            try:
                if not len(step) == len(self.distances):
                    raise ValueError("Inconsistent input for restraint generation: number of kappas in one or more timesteps is not the same as number of distances.")
            except:
                raise ValueError("Inconsistent input for restraint generation: kappas is not list of lists.")
        # Update internal values
        self.time_series = time_series
        self.kappas = kappas

    def add_distance_values(self, distance_series: list):
        """Manually add a series of distance values to be used by the PLUMED moving restraint.

        :param distance_series:     List of lists to record all the target value for the restraint. First dimension are time-steps, second dimension are distance values.
        :type distance_series:      List<List<float>>
        :raises ValueError:         Inconsistent input for restraint generation (some consistency checks)
        :raises RuntimeError:       Need to add timeseries before adding distance values.
        """
        # Perform consistency checks
        try:
            if not len(distance_series) == len(self.time_series):
                raise ValueError("Inconsistent input for restraint generation: lenghts of time_series and distance_series do not match.")
            for step in distance_series:
                try:
                    if not len(step) == len(self.distances):
                        raise ValueError("Inconsistent input for restraint generation: number of target distances in one or more timesteps is not the same as number of saved distances.")
                except:
                    raise ValueError("Inconsistent input for restraint generation: distance_series is not list of lists.")
        except AttributeError:
            raise RuntimeError("Need to add timeseries before adding distance values.")
        # Update internal values
        self.distance_series = distance_series

    def determine_distance_values(self, universe:MDAnalysis.Universe, frames: list, trj_time_step = 50):
        """Automatically add a series of distance values to be used by the PLUMED moving restraint, determined by averages over frames of the given trajectory.

        :param universe:        Universe which holds the trajectory.
        :type universe:         MDAnalysis.Universe
        :param frames:          List of 2-tuples, which holds the start- end points (in ps) to be used for averaging for each restraint step.
        :type frames:           List<Tuple<int>>
        :param trj_time_step:   Length of a frame of the trajectory in ps, defaults to 50.
        :type trj_time_step:    float, optional
        :raises RuntimeError:   Need to add timeseries before adding distance values.
        :raises ValueError:     Inconsistent input for restraint generation (some consistency checks)
        """
        # Consistency checks
        try:
            self.time_series
        except AttributeError:
            raise RuntimeError("Need to add timeseries before adding distance values.")
        if not len(frames) == len(self.time_series):
            raise ValueError("Inconsistent input for restraint generation: number of restraint steps in distances and frames does not match.")
        for frame in frames:
            if not isinstance(frame, tuple):
                raise ValueError("Inconsistent input for restraint generation: frames is not list of tuples.")
        # Determine distances from trajectory by averagign over the given frame ranges
        self.distance_series = []
        for frame in frames:
            startpoint = int(frame[0] / 50)
            endpoint = int(frame[1] / 50)
            distance_values = np.zeros((endpoint-startpoint, len(self.distances)))
            for i in range(0,endpoint-startpoint):
                # Shift the trajectory to the appropriate time point and find distances
                universe.trajectory[i+startpoint]
                distance_values[i,:] = np.array([dist.get_distance() for dist in self.distances])
            #Average those distances for this frame
            self.distance_series.append([np.mean(distance_values[:,i]) for i in range(len(self.distances))])
                


    def get_plumed_str(self, decimals = 4):
        """Obtain the PLUMED representation of this PLUMED_Restraint as a list of strings (multiline command).

        :raises RuntimeError:   Need to add timeseries and target values before generating the restraint string.
        :return:                List of strings (lines) to define the PLUMED restraint.
        :rtype:                 List<str>
        """
        # Check that all required fields are defined
        try:
            self.time_series
            self.distance_series
        except:
            raise RuntimeError("Need to add timeseries and target values before generating the restraint string.")
        # Create the list of strings which define the PLUMED moving restraint
        output = ["restraint: ...", "     MOVINGRESTRAINT"]
        output.append("     ARG=" + ','.join([d.label for d in self.distances]))
        for i in range(len(self.time_series)):
            at, step, kappa = f"     AT{i}=", f"STEP{i}=", f"KAPPA{i}="
            output.append(at + ','.join([str(np.round(x, decimals=decimals)) for x in self.distance_series[i]]) + "   " + \
                step + str(self.time_series[i]) + "   " + \
                kappa + ','.join([str(x) for x in self.kappas[i]]))
        output.append("...")
        return output