""" This holds the code for the main handler class PluMDAnalysis """

# Necessary imports
import MDAnalysis
from MDAnalysis.core.groups import AtomGroup

from .containers import PLUMED_Group
from .single_value_fields import PLUMED_Distance, PLUMED_Combine
from .plumed_restraint import PLUMED_Restraint

class PluMDAnalysis:
        """ Class to represent the PLUMED input file which is in production. Contains methods for adding elements and linking them,
        and will eventually output the PLUMED input file. 
        """
        def __init__(self, universe: MDAnalysis.Universe, trj_time_step = 50):
                """Initiate an instance of PluMDAnalysis, with an MDAnalysis Universe object holding the trajectory from which the PLUMED
                input file is to be generated. 

                :param universe:        MDAnalysis.Universe instance holding the trajectory on which the PLUMED input file will be based.
                :type universe:         MDAnalysis.Universe
                :param trj_time_step:   Length of a frame of the trajectory in ps, defaults to 50.
                :type trj_time_step:    float, optional
                :raises TypeError:      Constructor argument should be of type MDAnalysis.Universe.
                """        
                if not isinstance(universe, MDAnalysis.Universe):
                        raise TypeError("Constructor argument should be of type MDAnalysis.Universe.")
                self.plumed_groups = []
                self.single_value_fields = []
                self.restraints = []
                self.universe = universe
                self.trj_time_step = trj_time_step

        def add_atom_group(self, atom_group: AtomGroup, group_type = 'CENTER'):
                """Add an AtomGroup to this input script.

                :param atom_group:      AtomGroup to be added
                :type atom_group:       AtomGroup
                :param group_type:      Type of PLUMED group, 'CENTER' or 'COM', defaults to 'CENTER'
                :type group_type:       str, optional
                :return:                The created PLUMED_group
                :rtype:                 PLUMED_Group
                """
                plumed_group = self._contains_atom_group(atom_group)

                if not plumed_group:
                        plumed_group = PLUMED_Group(atom_group, f"group{len(self.plumed_groups) + 1}", group_type=group_type)
                        self.plumed_groups.append(plumed_group)
                return(plumed_group)

        def _contains_atom_group(self, atom_group: AtomGroup):
                """Determines whether an atom_group is already represented by a PLUMED_Group in this input script.

                :param atom_group:      AtomGroup to be tried. 
                :type atom_group:       AtomGroup
                :return:                Returns the corresponding PLUMED_Group if appropriate, or None if not contained.
                :rtype:                 PLUMED_Group or None
                """
                for plumed_group in self.plumed_groups:
                        if plumed_group.atom_group == atom_group:
                                return plumed_group
                return None
                

        def add_distance(self, atom_group_1: AtomGroup, atom_group_2: AtomGroup):
                """Add a new distance to this input script.

                :param atom_group_1:    AtomGroup for the first endpoint to the distance.
                :type atom_group_1:     AtomGroup
                :param atom_group_2:    AtomGroup for the second endpoint to the distance.
                :type atom_group_2:     AtomGroup
                :return:                The created PLUMED_Distance
                :rtype:                 PLUMED_Distance
                """
                # Check if the corresponding PLUMED_Group objects already exist. If yes, take them. Otherwise, create new ones.
                plumed_group_1 = self._contains_atom_group(atom_group_1)
                if not plumed_group_1:
                        plumed_group_1 = self.add_atom_group(atom_group_1)

                plumed_group_2 = self._contains_atom_group(atom_group_2)
                if not plumed_group_2:
                        plumed_group_2 = self.add_atom_group(atom_group_2)

                distance = PLUMED_Distance(plumed_group_1, plumed_group_2, f"distance{len(list(filter(lambda x: isinstance(x, PLUMED_Distance), self.single_value_fields))) + 1}")
                self.single_value_fields.append(distance)
                return distance

        def add_linear_combination(self, coefficients:list, fields_of_interest:list = ['all']):
                if fields_of_interest[0] == 'all':
                        new_combination = PLUMED_Combine(self.single_value_fields.copy(), coefficients, f"combine{len(list(filter(lambda x: isinstance(x, PLUMED_Combine), self.single_value_fields))) + 1}")
                else:
                        for field in fields_of_interest:
                                if not field in self.single_value_fields:
                                        raise ValueError("Provided fields_of_interest is not a subset of those already defined.")
                        new_combination = PLUMED_Combine(fields_of_interest, coefficients,f"combine{len(list(filter(lambda x: isinstance(x, PLUMED_Combine), self.single_value_fields))) + 1}")
                self.single_value_fields.append(new_combination)
                return new_combination

        def add_restraint_manual(self, time_series:list, kappas:list, series:list, fields_of_interest:list = ['all']):
                """Creates a new restraint from manual entry of the target single_value_fields at each step.

                :param time_series:         List of time-series points at which the PLUMED moving restraint takes steps.
                :type time_series:          List<int>
                :param kappas:              List of lists to record all the kappa values required for the restraint. First dimension are time-steps, second dimension are the single_value_fields.
                :type kappas:               List<List<int>>
                :param series:              List of lists to record all the target value for the restraint. First dimension are time-steps, second dimension are distance values.
                :type series:               List<List<float>>
                :param fields_of_interest:  List of all instances of SingleValueField to be considered for this restraint, subset of those already defined, ['all'] will use all which have been stored in the class, defaults to ['all']
                :type fields_of_interest:   List<SingleValueField>, optional
                :return:                    The created PLUMED_Restraint
                :rtype:                     PLUMED_Restraint
                """
                if fields_of_interest[0] == 'all':
                        new_restraint = PLUMED_Restraint(self.single_value_fields)
                else:
                        for field in fields_of_interest:
                                if not field in self.single_value_fields:
                                        raise ValueError("Provided fields_of_interest is not a subset of those already defined.")
                        new_restraint = PLUMED_Restraint(fields_of_interest)
                new_restraint.add_time_series(time_series, kappas)
                new_restraint.add_values(series)
                self.restraints.append(new_restraint)
                return new_restraint

        def add_restraint_from_trajectory(self, time_series:list, kappas:list, frames:list, fields_of_interest:list = ['all']):
                """Create a new restraint from automatically averaging all saved single_value_fields across the specified frames of the trajectory.

                :param time_series:         List of time-series points at which the PLUMED moving restraint takes steps.
                :type time_series:          List<int>
                :param kappas:              List of lists to record all the kappa values required for the restraint. First dimension are time-steps, second dimension are the single_value_fields.
                :type kappas:               List<List<int>>
                :param frames:              List of 2-tuples, which holds the start- end points (in ps) to be used for averaging for each restraint step.
                :type frames:               List<Tuple<int>>
                :param fields_of_interest:  List of all instances of SingleValueField to be considered for this restraint, subset of those already defined, ['all'] will use all which have been stored in the class, defaults to ['all']
                :type fields_of_interest:   List<SingleValueField>, optional
                :return:                    The created PLUMED_Restraint
                :rtype:                     PLUMED_Restraint
                """
                if fields_of_interest[0] == 'all':
                        new_restraint = PLUMED_Restraint(self.single_value_fields)
                else:
                        for field in fields_of_interest:
                                if not field in self.single_value_fields:
                                        raise ValueError("Provided fields_of_interest is not a subset of those already defined.")
                        new_restraint = PLUMED_Restraint(fields_of_interest)
                new_restraint.add_time_series(time_series, kappas)
                new_restraint.determine_values(self.universe,frames, trj_time_step=self.trj_time_step)
                self.restraints.append(new_restraint)
                return new_restraint

        def generate_PLUMED_input(self, print_args = [25000,'*','COLVAR'], write_to_file = True, filename = "plumed.dat"):
                """Compile the PLUMED input file out of everything added so far. Write to file if desired. 

                :param print_args:              Arguments for the PLUMED PRINT command, [STRIDE, ARG, FILE]. Defaults to [25000,'*','COLVAR']
                :type print_args:               list, optional
                :param write_to_file:           Whether or not the output should be written to file, defaults to True
                :type write_to_file:            bool, optional
                :param filename:                Filename for the input file to be written, defaults to "plumed.dat"
                :type filename:                 str, optional
                :return:                        The PLUMED input file as list of strings (lines).
                :rtype:                         List<str>
                """
                # build up strings from everything which was added
                output = []

                for group in self.plumed_groups:
                        output.append(group.get_plumed_str())
                output.append("")

                for distance in self.single_value_fields:
                        output.append(distance.get_plumed_str())
                output.append("")

                for restraint in self.restraints:
                        # remember that restraint gives a list, not a single string
                        output += restraint.get_plumed_str()
                output.append("")

                output.append(f"PRINT STRIDE={print_args[0]} ARG={print_args[1]} FILE={print_args[2]}")

                # Write to plumed output file if desired
                if write_to_file:
                        with open(filename, 'w') as f:
                                for line in output:
                                        f.write(line+'\n')

                return output
