from esi import calibration

dir = '../../Raw/'
prefix = 'e171118_'
hgne = '0004' 
xe = '0006' 
cuar = '0010'
flat = '0011,0012,0013' # dome flat
star = '0014,0015' # pinhole flat
bias = '0016,0017,0018,0019,0020,0021,0022,0023,0024,0025' # bias
out = 'calib'

calibration.prepare(dir,prefix,bias,star,hgne,cuar,xe,flat,out,onearc=False)
#calibration.prepare(dir,prefix,bias,star,hgne,cuar,xe,flat,out,arc=cuar,
#                    onearc=True)
# If arcs were taken separately, set onearc=False and do NOT use the arc
# parameter (which in the above example is set to arc=cuar)
