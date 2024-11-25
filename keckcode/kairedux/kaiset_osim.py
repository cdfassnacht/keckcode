"""

kaiset_osim.py

Define a class to run KAI functions on OSIRIS imaging data
"""

import numpy
from ..osiris.osimset import OsImSet

pyversion = sys.version_info.major

# ---------------------------------------------------------------------------


class KaiSetOsIm(OsImSet):
    """

    Class to run KAI functions within an OsImSet structure

    """

    def __init__(self, inlist, wcstype='koa', **kwargs):

        """ Make sure that inlist is in the correct format """
        if isinstance(inlist, (list, tuple, dict)):
            pass
        else:
            raise TypeError('\nKaiSetOsIm: inlist must be either a list, a'
                            ' tuple, or a dict')

        """ Call the superclass """
        if pyversion == 2:
            super(KaiSetOsIm, self).__init__(inlist, wcstype=wcstype, **kwargs)
        else:
            super().__init__(inlist, wcstype=wcstype, **kwargs)


