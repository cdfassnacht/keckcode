from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0739+1350'
name='S0739'
frames = [55,56]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 1.1 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
