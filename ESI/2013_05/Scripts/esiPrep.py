from esi import calibration

dir = '../Raw/'
prefix = 'e130514_'
cuar = '2001'
hgne = '1237' 
xe = '1239' 
flat = '0006,0007,0008,0009,0010,0011,0012,0013,0014,0015' # dome flat
star = '0016,0017,0018' # pinhole flat
bias = '0019,0020,0021,0022,0023,0024,0025,0026,0027,0028' # bias
out = 'calib'

calibration.prepare(dir,prefix,bias,star,hgne,cuar,xe,flat,out,arc=cuar,
                    onearc=True)
# If arcs were taken separately, set onearc=False and do NOT use the arc
# parameter (which in the above example is set to arc=cuar)
