from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S1010+5705'
name='S1010'
frames = [65,66]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 1.3 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
    #ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
