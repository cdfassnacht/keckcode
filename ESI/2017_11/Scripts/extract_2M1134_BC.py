from extract_generic import extract

#method = 'cdf'
method = 'oldham'
stdOrderCorr = 'orderCorr_Feige 110.dat'
fullname = '2M1134-2103BC'
name='2M1134_BC'
frames = [61, 62, 63]
#apcent = [0.,]
#aplab = ['test']
apcent = [-1.67, -0.35, 0, 2.0]
aplab = ['B', 'G', 'A', 'C'] # Check order of B and C
nsig = 0.8 # Width of aperture in terms of sigma.  Normal value is 1.0
apmin = -1.*nsig
apmax = nsig

for i in range(len(aplab)):
    extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig,
            method=method,apmin=apmin,apmax=apmax)
