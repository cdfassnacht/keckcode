from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S2211+1929'
name='S2211'
frames = [42,43,44]
apcent = [0.3,-0.69]
aplab = ['A','B']
nsig = 0.7 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
