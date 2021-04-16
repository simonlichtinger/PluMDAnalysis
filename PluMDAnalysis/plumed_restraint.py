""" This holds the PLUMED_Restraint class to represent a moving restraint. """

import MDAnalysis
import numpy as np

class PLUMED_Restraint:
    """ Class to define a PLUMED moving restraint, and output it as a list of strings. """
    def __init__(self, single_value_fields: list):
        """Constructor with the SingleValueField objects to be included in the restraint.

        :param single_value_fields:   List of SingleValueField objects to be included in the restraint.
        :type single_value_fields:    List<SingleValueField>
        """
        self.single_value_fields = single_value_fields

    def add_time_series(self, time_series: list, kappas: list):
        """Adds data for time-steps and kappa values across which the restraint is to work.

        :param time_series:  List of time-series points at which the PLUMED moving restraint takes steps.
        :type time_series:   List<int>
        :param kappas:       List of lists to record all the kappa values required for the restraint. First dimension are time-steps, second dimension are the single_value_fields.
        :type kappas:        List<List<int>>
        :raises ValueError:  Inconsistent input for restraint generation (some consistency checks)
        """
        # Do some consistency checks
        if not len(time_series) == len(kappas):
            raise ValueError("Inconsistent input for restraint generation: lenghts of time_series and kappa do not match.")
        for step in kappas:
            try:
                if not len(step) == len(self.single_value_fields):
                    raise ValueError("Inconsistent input for restraint generation: number of kappas in one or more timesteps is not the same as number of single_value_fields.")
            except:
                raise ValueError("Inconsistent input for restraint generation: kappas is not list of lists.")
        # Update internal values
        self.time_series = time_series
        self.kappas = kappas

    def add_values(self, series: list):
        """Manually add a series of distance values to be used by the PLUMED moving restraint.

        :param series:     List of lists to record all the target value for the restraint. First dimension are time-steps, second dimension are distance values.
        :type series:      List<List<float>>
        :raises ValueError:         Inconsistent input for restraint generation (some consistency checks)
        :raises RuntimeError:       Need to add timeseries before adding distance values.
        """
        # Perform consistency checks
        try:
            if not len(series) == len(self.time_series):
                raise ValueError("Inconsistent input for restraint generation: lenghts of time_series and series do not match.")
            for step in series:
                try:
                    if not len(step) == len(self.single_value_fields):
                        raise ValueError("Inconsistent input for restraint generation: number of target single_value_fields in one or more timesteps is not the same as number of saved single_value_fields.")
                except:
                    raise ValueError("Inconsistent input for restraint generation: series is not list of lists.")
        except AttributeError:
            raise RuntimeError("Need to add timeseries before adding distance values.")
        # Update internal values
        self.series = series

    def determine_values(self, universe:MDAnalysis.Universe, frames: list, trj_time_step = 50):
        """Automatically add a series ofvalues to be used by the PLUMED moving restraint, determined by averages over frames of the given trajectory. 

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
            raise ValueError("Inconsistent input for restraint generation: number of restraint steps in single_value_fields and frames does not match.")
        for frame in frames:
            if not isinstance(frame, tuple):
                raise ValueError("Inconsistent input for restraint generation: frames is not list of tuples.")
        # Determine single_value_fields from trajectory by averagign over the given frame ranges
        self.series = []
        for frame in frames:
            startpoint = int(frame[0] / 50)
            endpoint = int(frame[1] / 50)
            values = np.zeros((endpoint-startpoint, len(self.single_value_fields)))
            for i in range(0,endpoint-startpoint):
                # Shift the trajectory to the appropriate time point and find single_value_fields
                universe.trajectory[i+startpoint]
                values[i,:] = np.array([dist.get_value() for dist in self.single_value_fields])
            #Average those single_value_fields for this frame
            self.series.append([np.mean(values[:,i]) for i in range(len(self.single_value_fields))])
                


    def get_plumed_str(self, decimals = 4):
        """Obtain the PLUMED representation of this PLUMED_Restraint as a list of strings (multiline command).

        :raises RuntimeError:   Need to add timeseries and target values before generating the restraint string.
        :return:                List of strings (lines) to define the PLUMED restraint.
        :rtype:                 List<str>
        """
        # Check that all required fields are defined
        try:
            self.time_series
            self.series
        except:
            raise RuntimeError("Need to add timeseries and target values before generating the restraint string.")
        # Create the list of strings which define the PLUMED moving restraint
        output = ["restraint: ...", "     MOVINGRESTRAINT"]
        output.append("     ARG=" + ','.join([d.label for d in self.single_value_fields]))
        for i in range(len(self.time_series)):
            at, step, kappa = f"     AT{i}=", f"STEP{i}=", f"KAPPA{i}="
            output.append(at + ','.join([str(np.round(x, decimals=decimals)) for x in self.series[i]]) + "   " + \
                step + str(self.time_series[i]) + "   " + \
                kappa + ','.join([str(x) for x in self.kappas[i]]))
        output.append("...")
        return output