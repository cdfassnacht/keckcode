from esi import calibration

dir = '../Raw/'
prefix = 'e160708_'
cuar = '0015'
hgne = '0009' 
xe = '0012' 
flat = '0018,0019,0020' # dome flat
star = '0021,0022' # pinhole flat
bias = '0023,0024,0025,0026,0027,0028,0029,0030,0031,0032' # bias
out = 'calib'

calibration.prepare(dir,prefix,bias,star,hgne,cuar,xe,flat,out,onearc=False)
#calibration.prepare(dir,prefix,bias,star,hgne,cuar,xe,flat,out,arc=cuar,
#                    onearc=True)
# If arcs were taken separately, set onearc=False and do NOT use the arc
# parameter (which in the above example is set to arc=cuar)
