from lris.lris_pipeline import lris_pipeline

rawdir = '../Raw/'
prefix = 'lred'
sciframes  = '0121,0124,0125,0128,0129'
arcframe   = '0123'
flatframes = '0122,0126'
use_arc = 0
use_flat = 0
cache = 0
outpref = '0435m1c_red'

lris_pipeline(prefix,rawdir,sciframes,arcframe,flatframes,outpref,
              use_flat,use_arc,cache)
