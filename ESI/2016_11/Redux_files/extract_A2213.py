from extract_generic import extract

stdOrderCorr = 'orderCorr_Feige110.dat'
fullname = 'A2213-2652'
name='A2213'
frames = [39,40,41]
apcent = [0.,]
aplab = ['all',]

for jj in range(len(aplab)):
    ow,s,v = extract(fullname,name,frames,jj,apcent,aplab,stdOrderCorr,wid=1.0)
