"""
Script to make plots of final F0001 spectra
"""

from matplotlib import pyplot as plt
import spec_simple as ss

""" 
Read in the reduced spectra and add an offset to the flux of the brighter one
"""
F0001_A = ss.Spec1d('F0001_spec_ap_A.fits',informat='fits',logwav=True)
F0001_B = ss.Spec1d('F0001_spec_ap_B.fits',informat='fits',logwav=True)
tmpflux = F0001_A.flux + 0.5

""" Smooth the B spectrum """
F0001_B.smooth_boxcar(11,doplot=False)

""" Plot the spectra of the two components"""
plt.figure(1)
plt.plot(F0001_A.wav,tmpflux,'b')
plt.plot(F0001_B.wav,F0001_B.smoflux,'g')
plt.xlim(4400.,9250.)
plt.ylim(-0.15,2.)
z = 7167./1909. - 1.
ss.mark_spec_emission(z)
plt.title('F0001-0735 ESI Spectra 2016_07')

""" Mark the lensing galaxy spectral lines """
plt.figure(2)
plt.plot(F0001_B.wav,F0001_B.smoflux,'g')
plt.xlim(4400.,9250.)
plt.ylim(-0.15,0.4)
F0001_B.lineinfo['dxlab'][F0001_B.lineinfo['name']=='CaII H'] = 40.
F0001_B.lineinfo['dxlab'][F0001_B.lineinfo['name']=='G-band'] = -20.
F0001_B.lineinfo['dxlab'][F0001_B.lineinfo['name']=='H-gamma'] = 10.
F0001_B.lineinfo['dxlab'][F0001_B.lineinfo['name']=='Fe4383'] = 35.
F0001_B.mark_spec_absorption(0.69,usesmooth=True,labloc='topright',tickfac=0.4)
plt.title('F0001-0735 Lensing Galaxy Lines')
plt.show()
