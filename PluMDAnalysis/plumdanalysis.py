""" This holds the code for the main handler class PluMDAnalysis """

# Necessary imports
import MDAnalysis as mda

class PluMDAnalysis:
        """ Class to represent the PLUMED input file which is in production. Contains methods for adding elements and linking them,
        and will eventually output the PLUMED input file. 
        """
        def __init__(self, universe: mda.Universe):
                """Initiate an instance of PluMDAnalysis, with an MDAnalysis Universe object holding the trajectory from which the PLUMED
                input file is to be generated. 

                :param universe: MDAnalysis.Universe instance holding the trajectory on which the PLUMED input file will be based.
                :type universe: MDAnalysis.Universe
                :raises TypeError: Constructor argument should be of type MDAnalysis.Universe.
                """        
                if not type(universe) == mda.Universe:
                        raise TypeError("Constructor argument should be of type MDAnalysis.Universe.")
                self.universe = universe
        




