from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0011-0845'
name='S0011'
frames = [35,36]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 2.5 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
