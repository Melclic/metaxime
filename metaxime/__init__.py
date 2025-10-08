from .version import __version__

#from . import utils
#from .parser import ParserRP2
#from .cache_data import RR_Data

#make sure that RDKit does not have those warning messages
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)   # only show critical messages

import logging

logging.basicConfig(
    #level=logging.DEBUG,
    #level=logging.INFO,
    #level=logging.WARNING,
    level=logging.ERROR,
    format='%(asctime)s.%(msecs)03d %(levelname)s %(module)s - %(funcName)s: %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S',
)

__all__ = ["__version__"]
