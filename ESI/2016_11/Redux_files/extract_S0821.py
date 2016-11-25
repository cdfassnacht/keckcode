from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0821+3404'
name='S0821'
frames = [59,60]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 1.4 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
