from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0748+4225'
name='S0748'
frames = [56,57]
apcent = [0.,]
aplab = ['all',]
nsig = 1.2 # Width of aperture in terms of sigma.  Normal value is 1.0

for i in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig)
