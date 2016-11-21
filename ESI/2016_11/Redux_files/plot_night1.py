"""
Plots the total light spectra from night 1 in one figure
"""

from astropy.io import ascii
import spec_simple as ss
from matplotlib import pyplot as plt


lensinfo = ascii.read('redshifts.txt')
plt.figure(figsize=(6,9))
plt.subplots_adjust(hspace=0.001)
count = 1
axlist = []
for i in lensinfo:
    infile = '%s_spec_ap_all.fits' % i['Lens']
    print infile
    spec = ss.Spec1d(infile,informat='fits',logwav=True)
    if count == 1:
        ax = plt.subplot(8,1,count)
    else:
        ax = plt.subplot(8,1,count,sharex=axlist[0])
    axlist.append(ax)
    if i['smo'] > 1:
        spec.smooth_boxcar(i['smo'],title=None)
        smo = True
    else:
        spec.plot(title=None)
        smo = False
    plt.xlim(4150.,10300.)
    plt.ylim(-0.03,i['ymax'])
    spec.mark_speclines('em',i['zs'],usesmooth=smo,showz=False)
    del(spec)
    count += 1

xticklabels = axlist[0].get_xticklabels()
yticklabels = axlist[0].get_yticklabels()
for i in range(1,7):
    xticklabels += axlist[i].get_xticklabels()
    yticklabels += axlist[i].get_yticklabels()
plt.setp(xticklabels, visible=False)
plt.setp(yticklabels, visible=False)

plt.show()
