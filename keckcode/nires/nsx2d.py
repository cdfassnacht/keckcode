"""

nsx2d.py

A class written to take the *corrected.fits files produced by running
Tom Barlow's NSX program to reduce NIRES data, and to process them as
echelle2d.Ech2d structures.

"""

import numpy as np
from astropy.io import fits as pf
from specim import imfuncs as imf
from specim.specfuncs import echelle2d as ech2d


# ---------------------------------------------------------------------------

class Nsx2d(ech2d.Ech2d):
    """

    A class written to take the *corrected.fits files produced by running
    Tom Barlow's NSX program to reduce NIRES data, and to process them as
    echelle2d.Ech2d structures.

    """

    def __init__(self, inroot, frame, frame2=None, scisuff='corrected',
                 varsuff=None):
        """

        The class gets instantiated by reading in the single fits file that
        contains the five NIRES order in one HDU image, and the associated
        sky file, and makes five postage-stamp cutouts that can the be
        converted into ImageHDUs and then a single Ech2d object.

        """

        """ Set the parameters of the postage stamp regions """
        dtype = [('xmin', int), ('xmax', int), ('ymin', int), ('ymax', int)]
        oinfo = np.array([
                (13, 1005, 0, 115),
                (13, 2034, 200, 315),
                (13, 2034, 400, 515),
                (13, 2034, 600, 715),
                (13, 2034, 800, 915),
                ], dtype=dtype)

        """ Set the input filename(s) """
        if frame2 is not None:
            fileroot = '%s_%04d-%04d' % (inroot, frame, frame2)
        else:
            fileroot = '%s_%04d' % (inroot, frame)
        infile = '%s-%s.fits' % (fileroot, scisuff)
        if varsuff is not None:
            varfile = '%s-%s.fits' % (fileroot, varsuff)
        else:
            varfile = None
            varhdu = None

        """ Read the input science file into an Image class """
        print(infile)
        sci = imf.Image(infile)

        """ Make the cutouts """
        scihdu = pf.HDUList(pf.PrimaryHDU())
        vardata = []
        for info in oinfo:
            sci.set_subim_xy(info['xmin'], info['ymin'], info['xmax'],
                             info['ymax'], verbose=False)
            scihdu.append(pf.ImageHDU(sci.data, sci.subimhdr))

        """ Now put the data into the superclass """
        super(Nsx2d, self).__init__(scihdu, varhdu)

        """ Set the default multi-order grid """
        self.plotgrid = (1, 5)
