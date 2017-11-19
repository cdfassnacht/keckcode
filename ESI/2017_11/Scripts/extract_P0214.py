from extract_generic import extract

#method = 'cdf'
method = 'oldham'
stdOrderCorr = 'orderCorr_Feige 110.dat'
fullname = 'P0214+4618'
name='P0214'
frames = [42, 43]
#apcent = [0.,]
#aplab = ['all']
nsig = 1. # Width of aperture in terms of sigma.  Normal value is 1.0
apcent = [0, 2.16]
aplab = ['A', 'B']
nsig = 1. # Width of aperture in terms of sigma.  Normal value is 1.0
apmin = -1.*nsig
apmax = nsig

for i in range(len(aplab)):
    extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig,
            method=method,apmin=apmin,apmax=apmax)
