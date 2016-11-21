from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'A2338-2700'
name='A2338'
frames = [33,34]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 1. # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
