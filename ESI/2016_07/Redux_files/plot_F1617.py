"""
Script to make plots of final F1617 spectra
"""

from matplotlib import pyplot as plt
import spec_simple as ss

""" 
Read in the reduced spectra and add an offset to the flux of the brighter one
"""
F1617_A = ss.Spec1d('F1617_spec_ap_A.fits',informat='fits',logwav=True)
F1617_B = ss.Spec1d('F1617_spec_ap_B.fits',informat='fits',logwav=True)
tmpflux = F1617_A.flux + 0.5

""" Smooth the B spectrum """
F1617_B.smooth_boxcar(11,doplot=False)

""" Plot the spectra of the two components"""
plt.figure(1)
plt.plot(F1617_A.wav,tmpflux,'b')
plt.plot(F1617_B.wav,F1617_B.smoflux,'g')
plt.xlim(4400.,9250.)
plt.ylim(-0.1,1.)
z = 4769./1549. - 1.
ss.mark_spec_emission(z)
plt.title('F1617+3827 ESI Spectra 2016_07')

""" Mark the lensing galaxy spectral lines """
#plt.figure(2)
#plt.plot(F1617_B.wav,F1617_B.smoflux,'g')
#plt.xlim(4400.,9250.)
#plt.ylim(-0.15,0.4)
#F1617_B.lineinfo['dxlab'][F1617_B.lineinfo['name']=='CaII H'] = 40.
#F1617_B.lineinfo['dxlab'][F1617_B.lineinfo['name']=='G-band'] = -20.
#F1617_B.lineinfo['dxlab'][F1617_B.lineinfo['name']=='H-gamma'] = 10.
#F1617_B.lineinfo['dxlab'][F1617_B.lineinfo['name']=='Fe4383'] = 35.
#F1617_B.mark_spec_absorption(0.69,usesmooth=True,labloc='topright',tickfac=0.4)

plt.show()
