from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0940-0249'
name='S0940'
frames = [58,59]
apcent = [0.,-2.13]
#aplab = ['all']
aplab = ['A','B']
nsig = 0.7 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
