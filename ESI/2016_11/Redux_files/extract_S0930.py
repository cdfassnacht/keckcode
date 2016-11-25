from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0930+4614'
name='S0930'
frames = [65,]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 0.9 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
