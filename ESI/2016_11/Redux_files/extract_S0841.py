from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0841+0312'
name='S0841'
frames = [61,62]
apcent = [0.,2.96]
#aplab = ['all']
aplab = ['A','B']
nsig = 0.9 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
    #ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
