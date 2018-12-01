import pylab,scipy
import special_functions as sf
import matplotlib as mpl

class SpecID:

	def __init__(self,data):
		self.data = data.astype(scipy.float32)
		self.points = scipy.arange(data.size,dtype=scipy.float32)
		self.wave = self.points.copy()
		self.solution = None
		self.inverse = None
		self.ax = None
		self.lines = []
		self.markers = {}
		self.ids = []
		self.__call__()
		self.prompt = None
		self.pause = None

	def __call__(self):
		pylab.plot(self.wave,self.data)
		self.ax = pylab.gca()
		self.cid = pylab.connect('key_release_event',self.menu)

	def menu(self,event):
		if event.key=='m':
			self.fitline(event.xdata)
		elif event.key=='d':
			self.rmline(event.xdata)
		elif event.key=='s':
			self.dofit()

	def fitline(self,xpos):
		# Find local maximum to start fitting
		if self.solution==None:
			point = round(xpos)
		else:
			point = sf.genfunc(xpos,0.,self.inverse)
			point = point[0].round()
		print point
		center = self.data[point-5:point+6].argmax()+point-5
		print self.data[point-5:point+6]
		print self.data[point-5:point+6].argmax()

		max = self.data[center]
		fit = scipy.empty(4)
		fit[0] = 0.
		fit[1] = max
		fit[2] = 7.
		fit[3] = 2.

		fit,chi = sf.ngaussfit(self.data[center-7:center+8],fit)
		centroid = fit[2]+center-7
		print "Press enter to input wavelength"
		wave = float(input("Wavelength: "))
		max = self.data[centroid-2:centroid+3].max()
		try:
			indx = self.lines.index(centroid)
			self.ids[indx] = wave
		except:
			self.lines.append(centroid)
			self.ids.append(wave)
			axis = self.ax.axis()
			fudge = 0.05*(axis[3]-axis[2])
			if self.solution==None:
				wave = centroid
			else:
				wave = sf.genfunc(centroid,0.,self.solution)
			mark = mpl.patches.Polygon([(wave,max+fudge),(wave,max+2*fudge)])
			self.ax.add_patch(mark)
			self.markers[centroid] = mark
			pylab.draw()

	def rmline(self,xpos):
		if self.solution==None:
			point = xpos
		else:
			point = sf.genfunc(xpos,0.,self.inverse)
		delta = 1e7
		indx = -1
		for i in range(len(self.lines)):
			diff = abs(self.lines[i]-point)
			if diff<delta:
				delta = diff
				indx = i
		if indx>-1:
			val = self.lines[indx]
			self.lines.remove(val)
			self.ax.patches.remove(self.markers[val])
			del self.markers[val]
			val = self.ids[indx]
			self.ids.remove(val)
			pylab.draw()

	def dofit(self):
		print "Press Enter to input order"
		order = raw_input("Order: ")
		order = int(order)
		lines = scipy.empty(len(self.lines))
		wave = scipy.empty(len(self.ids))

		for i in range(lines.size):
			lines[i] = self.lines[i]
			wave[i] = self.ids[i]

		fit = scipy.empty((lines.size,2))
		fit[:,0] = lines.copy()
		fit[:,1] = wave.copy()
		self.solution = sf.lsqfit(fit,'chebyshev',order)
		fit[:,0] = wave.copy()
		fit[:,1] = lines.copy()
		self.inverse = sf.lsqfit(fit,'chebyshev',order)

		self.wave = sf.genfunc(self.points,0.,self.solution)
		pylab.close()
		self.__call__()
		for i in self.markers:
			wave = sf.genfunc(i,0.,self.solution)
			verts = self.markers[i].get_verts()
			bottom = verts[0][1]
			top = verts[1][1]
			wave = wave[0]
			mark = mpl.patches.Polygon([(wave,bottom),(wave,top)])
			self.markers[wave] = mark
			self.ax.add_patch(mark)
			del self.markers[i]
		pylab.draw()

