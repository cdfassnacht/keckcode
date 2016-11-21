from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'D0509-2350'
name='D0509'
frames = [52,53]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 2. # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
