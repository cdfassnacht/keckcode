from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0941+0518'
name='S0941'
frames = [60,61]
apcent = [0.,3.516,5.096]
aplab = ['A','G','B']
nsig = 1.0 # Width of aperture in terms of sigma.  Normal value is 1.0

for i in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig)
