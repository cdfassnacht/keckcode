import numpy as np
from astropy.table import Table
from specim.specfuncs import echelle1d

"""
============================== Esi1d class ==============================
"""


class Esi1d(echelle1d.Ech1d):
    """
    A class for ESI 1D spectra, which have been extracted by the Esi2d
    methods, but have not yet been combined into one final output spectrum.
    Therefore, there are 10 extracted 1d spectra, one for each order.
    These 10 extracted spectra will be stored in an array of Spec1d instances.

    """

    def __init__(self, inspec, informat='text', summary=True, verbose=True):
        """

        Initializes an Esi1d instance, essentially by initializing an
        Ech1d instance with the ESI order information

        """

        """
        Define the information pertaining to the ESI echelle orders

        Note that the pixmin and pixmax are not used at this point since
        any trimming of the orders should have been done in previous
        steps.
        """
        dtype = [('order', int), ('pixmin', int), ('pixmax', int)]
        oinfo = np.array([
                (1,  0, -1),
                (2,  0, -1),
                (3,  0, -1),
                (4,  0, -1),
                (5,  0, -1),
                (6,  0, -1),
                (7,  0, -1),
                (8,  0, -1),
                (9,  0, -1),
                (10, 0, -1),
                ], dtype=dtype)
        ordinfo = Table(oinfo)
        # ordinfo = None

        """ Initialize by calling the parent class """
        super(Esi1d, self).__init__(inspec, informat=informat, ordinfo=ordinfo,
                                    summary=summary, verbose=verbose)
