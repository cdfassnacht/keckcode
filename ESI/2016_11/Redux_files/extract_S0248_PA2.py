from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0248_PA2'
name='S0248'
frames = [41,42]
apcent = [0.,-1.18,3.9]
#aplab = ['all']
aplab = ['C','D','G2']
nsig = 0.8 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
