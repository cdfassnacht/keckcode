from esi import calibration

dir = '../Raw/'
prefix = 'e161120_'
cuar = '0019'
hgne = '0013' 
xe = '0015' 
flat = '0021,0022,0023' # dome flat
star = '0024,0025' # pinhole flat
bias = '0026,0027,0028,0029,0030,0031,0032,0033,0034,0035' # bias
out = 'calib'

calibration.prepare(dir,prefix,bias,star,hgne,cuar,xe,flat,out,onearc=False)
#calibration.prepare(dir,prefix,bias,star,hgne,cuar,xe,flat,out,arc=cuar,
#                    onearc=True)
# If arcs were taken separately, set onearc=False and do NOT use the arc
# parameter (which in the above example is set to arc=cuar)
