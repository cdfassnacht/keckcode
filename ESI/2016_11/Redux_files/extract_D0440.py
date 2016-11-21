from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'D0440-2008'
name='D0440'
frames = [50,51]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 2. # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
