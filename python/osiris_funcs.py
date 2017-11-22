"""

osiris_funcs.py

"""

import numpy as np
import imfuncs as imf

class osiris(imf.Image):
    """
    A class used to visualize and analyze OSIRIS data
    """

    def __init__(self, infile):
        """
        Loads in an OSIRIS data cube that was processed through the
        standard data reduction pipeline and, possibly, has had some additional
        processing done.
        """
