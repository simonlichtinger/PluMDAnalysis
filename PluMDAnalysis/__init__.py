""" This library automates the generation of PLUMED input files by providing an
interface with python trajectory analysis via MDAnalysis. The current focus is
on steered MD via distances between atom groups, but I will try to keep it 
extensible and will do so once the occasion arises. """

# Import version info
from .version_info import VERSION_INT, VERSION  # noqa

# Import PLUMED-handler class

from .plumdanalysis import PluMDAnalysis
