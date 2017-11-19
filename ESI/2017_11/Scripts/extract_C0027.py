from extract_generic import extract

#method = 'cdf'
method = 'oldham'
stdOrderCorr = 'orderCorr_Feige 110.dat'
fullname = 'C0027+0232'
name='C0027'
frames = [37, 38]
#apcent = [0.,]
#aplab = ['all']
#nsig = 1. # Width of aperture in terms of sigma.  Normal value is 1.0
apcent = [-2.73, 0]
aplab = ['B', 'A']
nsig = 1. # Width of aperture in terms of sigma.  Normal value is 1.0
apmin = -1.*nsig
apmax = nsig

for i in range(len(aplab)):
    extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig,
            method=method,apmin=apmin,apmax=apmax)
