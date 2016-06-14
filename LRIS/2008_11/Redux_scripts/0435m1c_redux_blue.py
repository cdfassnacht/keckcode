from lris.lris_pipeline import lris_pipeline

rawdir = '../Raw/'
prefix = 'lblue'
sciframes  = '0146,0151,0152,0156,0157'
arcframes  = '0154'
flatframes = '0147,0153,0158'
use_arc = 0
use_flat = 0
cache = 0
outpref = '0435m1c_blue'

lris_pipeline(prefix,rawdir,sciframes,arcframe,flatframes,outpref,
              use_flat,use_arc,cache)
