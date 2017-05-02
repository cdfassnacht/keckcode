from extract_generic import extract

stdOrderCorr = 'orderCorr_HZ44.dat'
fullname = 'EEL_J1605+3811'
name='J1605'
frames = [49, 50]
apcent = [0.,]
aplab = ['all']
# apcent = [0.,3.8]
# aplab = ['A','B']
nsig = 2. # Width of aperture in terms of sigma.  Normal value is 1.0

for i in range(len(aplab)):
    extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig)
