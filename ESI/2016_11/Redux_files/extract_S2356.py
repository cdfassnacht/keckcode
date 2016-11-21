from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S2356-0534'
name='S2356'
frames = [47,48]
apcent = [0.,]
aplab = ['all',]
nsig = 2.0 # Width of aperture in terms of sigma.  Normal value is 1.0

for i in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig)
