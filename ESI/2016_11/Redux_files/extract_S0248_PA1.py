from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0248_PA1'
name='S0248'
frames = [40,]
apcent = [0.,-1.3]
#aplab = ['all',]
aplab = ['A','B']
nsig = 1. # Width of aperture in terms of sigma.  Normal value is 1.0

for i in range(len(aplab)):
    extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig)
