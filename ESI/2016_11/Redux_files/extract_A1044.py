from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'A1044-1639'
name='A1044'
frames = [66,67]
apcent = [0.,]
aplab = ['all']
#aplab = ['A','B']
nsig = 1.7 # Width of aperture in terms of sigma.  Normal value is 1.0

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=nsig)
