from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'S0808-0051'
name='S0808'
frames = [54,55]
apcent = [0.,3.8]
aplab = ['A','B']
nsig = 1.2 # Width of aperture in terms of sigma.  Normal value is 1.0

for i in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,i,apcent,aplab,stdOrderCorr,wid=nsig)
