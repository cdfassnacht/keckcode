"""
LRIS reduction pipeline

The pipeline is generally accessed via a calling script that looks something
  like this:

	from lris import lris_pipeline

	dir = "data_directory/rawdata/"
	prefix = "lred" (or "lblue" or--for old data--"lris")
	outprefix = "mask2_11oct04"
	science = "0101,0102,0107,0108"
	arc = "0106"
	flats = "0103,0104,0105"
	usearc = 0
	useflat = 0
	cache = 0
	offsets = None

	lris_pipeline(prefix,dir,science,arc,flats,usearc,useflat,cache)

"""

import lris_blue,lris_red
import lris_pipeline,lris_biastrim
