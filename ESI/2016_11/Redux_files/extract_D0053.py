from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'D0053+0020'
name='D0053'
frames = [49,]
apcent = [0.,]
aplab = ['all',]
nsig = 1.8 # Width of aperture in terms of sigma.  Normal value is 1.0

for i in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig)
