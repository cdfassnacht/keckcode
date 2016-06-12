"""
Helper routines to perform bias subtraction and overscan trimming of LRIS data.
"""

import scipy


def oneamp(data):
	"""
	Subtracts bias from data and returns the overscan region-subtracted
	  image.
	"""
	bias = (data[:,2:21].mean(axis=1)*18+data[:,2069:2148].mean(axis=1)*80)/98.
	out_data = data[:,21:2069]-bias

	return out_data


def redside(data):
	"""
	Subtracts bias from data and returns the overscan region-subtracted
	  image. CCD geometry is currently hardwired, so this won't work for
	  windowed or binned setups.
	"""
	if data.shape[1]==2148:
		return oneamp(data)

	bias = scipy.empty((2,data.shape[0]))

	bias[0] = (data[:,1:20].mean(axis=1)*19+data[:,2088:2168].mean(axis=1)*80)/99.
	bias[1] = (data[:,20:38].mean(axis=1)*18+data[:,2168:2248].mean(axis=1)*80)/98.


	out_data = scipy.empty((2046,2048))

	"""
	Mask out the bad columns. Note this might not be appropriate for older
	  data (or if the CCDs change).
	"""
	mask = (data[0:1490,995]+data[0:1490,997])/2.
	data[0:1490,996] = mask.copy()
	mask = (data[0:1493,999]+data[0:1493,1001])/2.
	data[0:1493,1000] = mask.copy()

	data = data.transpose()

	out_data[0:1023,:] = data[41:1064,:] - bias[0]
	out_data[1023:2046,:] = data[1064:2087,:] - bias[1]

	"""
	Fix difference in amplifier gains. This does *not* convert from DN
	  to electrons.
	"""
	out_data[1023:2046,:] *= 1.0765 # Note this differs from the LRIS
					#  website that would suggest 1.0960

	out_data = out_data.transpose()

	return out_data

def blueside(data):
	"""
	Subtracts bias from data and returns the overscan region-subtracted
	  image.
	"""
	data = data.T
	data = data[:,::-1]
	size = data.shape[1]

	bias = scipy.empty((4,size))
	bias[0] = (data[0:50,:].mean(axis=0) + data[4300:4380,:].mean(axis=0))/2.
	bias[1] = (data[52:101,:].mean(axis=0) + data[4380:4460,:].mean(axis=0))/2.
	bias[2] = (data[102:153,:].mean(axis=0) + data[4460:4540,:].mean(axis=0))/2.
	bias[3] = (data[153:202,:].mean(axis=0) + data[4540:4620,:].mean(axis=0))/2.

	"""
	Conversion factor from DN to electrons (from LRIS website; disabled for
	  now).
	"""
	gain = [1.55,1.56,1.63,1.70]
	outdata = scipy.empty((4096,4096))
	for i in range(4):
		outstart = i*1024
		datastart = i*1024 + 204
		outdata[outstart:outstart+1024,:] = data[datastart:datastart+1024,:] - bias[i]
	#	outdata[outstart:outstart+1024,:] *= gain[i]
	del bias

	return outdata
