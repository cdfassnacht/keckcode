"""

nsxset.py - Code to perform actions on multiple nsxspec files

"""

from os import path
from matplotlib import pyplot as plt
from specim.specfuncs.specset1d import SpecSet1d
from .nsxspec import NsxSpec
from .nires1d import Nires1d


class NsxSet(list):
    """

    The NsxSet class is effectively just a list of NsxSpec instances that
    includes operations that act on the whole list (e.g., coadds)

    """

    # ------------------------------------------------------------------------

    def __init__(self, root, frames, frames2=None):

        """
        Reads the input data sets, which have been produced by running the
        NSX code on NIRES data.  The files are designated by a root string
        and frame numbers.  One of the options within NSX is to run the
        code simultaneously on both an A and B frame.  In that case, the
        extracted spectra will have a base name of root_frame_frame2

        Required inputs:
          root    - root name of files, e.g., 's190122'
          frames  - list of frame numbers (as integers)

        Optional inputs:
          frames2 - list of second frame numbers (as integers) if NSX was
                     run in A-B/B-A mode.  The default value (None) means that
                     NSX was NOT run in this mode
                    If set, then the base names will have the form
                     root_frame-frame2
                    If not set (default) then this list is not used and the
                     base names will have the form root_frame

        """

        """ Setup """
        if frames2 is None:
            frames2 = []
            for f in frames:
                f2 = None
                frames2.append(f2)

        """ Read in the data """
        for f, f2 in zip(frames, frames2):
            print('')
            d = NsxSpec(root, f, f2)
            self.append(d)

        """ Get the order list from the first input spectrum """
        self.ordinfo = self[0].ordinfo

    # ------------------------------------------------------------------------

    def coadd(self, doplot=True, outfile=None, smo=None, z=None, **kwargs):
        """

        Coadds the extracted 1d spectra from each order separately

        """

        """ Loop through the orders """
        coaddlist = []
        for order in range(len(self.ordinfo)):

            """ Create a list of Spec1d objects """
            speclist = []
            for nspec in self:
                if nspec[order] is not None:
                    speclist.append(nspec[order])
            if len(speclist) == 0:
                print('')
                print('ERROR: Called coadd but inputs do not have '
                      'extracted spectra yet')
                print('')
                raise ValueError

            """ Coadd the spectra in the list """
            specall = SpecSet1d(spec1dlist=speclist)
            coaddlist.append(specall.coadd(doplot=False, verbose=False,
                                           **kwargs))
            del(speclist)

        print('Reading coadded file')
        outspec = Nires1d(coaddlist)
        if doplot:
            outspec.plot_all(smo=smo, z=z, **kwargs)
        return outspec
