import optparse, re
import numpy as np
import math, scipy.stats
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


class CSM:
	def __init__(self):
		self.lookup = {}
		self.twoDslice = []
		self.minX = 0.0 
		self.minY = 0.0
		self.gridXsize = 0.0
		self.gridYsize = 0.0
		self.source_points = np.zeros(shape=(180,2))

	def projectScans(self, state, scans):
		scanAngles = np.array([minRange+increment*i for i in range(numScans)])
		#print scanAngles
		orientation = np.array([state[2],]*numScans)
		scan = np.array(scans, dtype=np.float32)
		#print scan.shape
		x = scan*np.cos(scanAngles+orientation)
		y = scan*np.sin(scanAngles+orientation)
		projectedPoints = np.array([state[0]+x, state[1]+y])
		return np.transpose(projectedPoints)

	def buildLookupTable(self, prev_state, prev_scans):
		mappoints = self.projectScans(prev_state, prev_scans)
		print "The shape of mappoints is "+str(mappoints.shape)
		maxX, maxY = tuple(mappoints.max(axis=0))
		minX, minY = tuple(mappoints.min(axis=0))
		maxX, maxY = math.ceil(maxX), math.ceil(maxY)
		minX, minY = math.floor(minX), math.floor(minY)
		gridRangeX = maxX-minX
		gridRangeY = maxY-minY
		#final map size by extrapolating the range to left, right, top and bottom
		#This is done as scan points might fall outside the original grid size
		self.minX = minX-gridRangeX
		self.maxX = maxX+gridRangeX
		self.minY = minY-gridRangeY
		self.maxY = maxY+gridRangeY
		gridRangeX = self.maxX-self.minX
		gridRangeY = self.maxY-self.minY
		numgrids = 1000
		self.gridXsize = float(gridRangeX)/numgrids
		self.gridYsize = float(gridRangeY)/numgrids
		lookup = {}
		for i in range(numgrids):
			for j in range(numgrids):
				x, y = self.minX+self.gridXsize*(i+0.5), self.minY+self.gridYsize*(j+0.5)
				self.lookup[(i, j)] = self.closestPoint(x, y, mappoints)

	def closestPoint(self, x, y, mappoints):
		dist = np.sqrt(np.square(mappoints[:, 0]-x)+np.square(mappoints[:, 1]-y))
		minDist = dist.min()
		return 1.0/minDist
		#return -minDist

	def computeGMean(self, source_points):
		weights = []
		for point in source_points:
			weights.append(self.lookup[tuple(point)])
		#mini = min(weights)
		#weights = np.array(weights)-mini
		#print max(weights)
		#prob = scipy.stats.mstats.gmean(weights)
		prob = sum(weights)/len(weights)
		#print prob
		return prob

	def twoDslices(self, scans, state):
		boxsize = 100
		slice = np.zeros(shape=(boxsize, boxsize))
		state = np.array(state)
		state.shape = (1,3)
		#Grid corresponding to the state
		gridCenter = self.findGrid(state)
		#Reset grid to the bottom leftmost position
		gridOriginCenter = gridCenter-[50*self.gridXsize, 50*self.gridYsize]
		gridCenter = gridOriginCenter
		prev_source_points = self.source_points
		self.source_points = self.projectScans(np.append(gridCenter, state[0][2]), scans)
		print self.source_points.shape
		self.source_points = self.findGrid(self.source_points)
		print self.source_points.shape
		if (prev_source_points==self.source_points).all():
			print "THEY ARE THE SAME FOR SOURCE POINTS"
		source_points_transformed = np.zeros(shape=self.source_points.shape)
		for i in range(100):
			for j in range(100):
				#source_points_transformed = source_points+[self.gridYsize*j, self.gridXsize*i]
				source_points_transformed = self.source_points+[j, i]
				slice[i][j] = self.computeGMean(source_points_transformed)
				#gridCenter = gridOriginCenter+[self.gridXsize*j, self.gridYsize*i]
		return slice

	def findGrid(self, state):
		index = np.array([list((state[:, 0]-self.minX)/self.gridXsize), list((state[:, 1]-self.minY)/self.gridYsize)])
		index = np.transpose(index)
		#print index.astype(int)
		#gridCenter = [self.minX, self.minY]+[self.gridXsize, self.gridYsize]*(index.astype(int)+0.5)
		return index.astype(int)


	def findGrid_(self, state):
		print state
		index = np.array([(state[0]-self.minX)/self.gridXsize, (state[1]-self.minY)/self.gridYsize])
		index = np.transpose(index)
		return index.astype(int)

def convert_to_radians(angle):
	return angle*math.pi/180

def plotGaussian(slice_num, state):
	slice = slices[slice_num]
	theta = -rad+(slice_num*2*rad/30)+orientation
	print state
	state = [state[0], state[1], theta]
	print state
	gridCenter = csm.findGrid_(state)
	gridOriginCenter = gridCenter-[50*csm.gridXsize, 50*csm.gridYsize]
	x_step = csm.gridXsize
	y_step = csm.gridYsize
	meanX = 0.0
	meanY = 0.0
	for i in range(100):
		x = x_step*i+gridOriginCenter[0]
		for j in range(100):
			meanX+=slice[i, j]*x
	for i in range(100):
		y = i*y_step+gridOriginCenter[1]
		for j in range(100):
			meanY+=slice[j, i]*y

	meanX = meanX/10000
	meanY = meanY/10000
	covrXX = 0.0
	covrXY = 0.0
	covrYY = 0.0
	for i in range(100):
		for j in range(100):
			x = gridOriginCenter[0]+x_step*i
			y = gridOriginCenter[1]+y_step*j
			x_ = gridOriginCenter[0]+x_step*j
			y_ = gridOriginCenter[1]+y_step*i
			covrXY+=(slice[i, j]*x-meanX)*(slice[i, j]*y-meanY)
			covrXX+=(slice[i, j]*x-meanX)*(slice[i, j]*x_-meanX)
			covrYY+=(slice[i, j]*y_-meanY)*(slice[i, j]*y-meanY)
	covrXY = covrXY/100000000
	covrXX = covrXX/100000000
	covrYY = covrYY/100000000
	print "MeanX is "+str(meanX)
	print "MeanY is "+str(meanY)
	print "covrXX is "+str(covrXX)
	print "covrXY is "+str(covrXY)
	print "covrYY is "+str(covrYY)
	NUM = 1
	ells = [Ellipse(xy=(11,37), width=covrXX, height=covrYY, angle=np.rad2deg(np.arccos(-0.45)))
	        for i in range(NUM)]

	fig = plt.figure(0)
	ax = fig.add_subplot(111, aspect='equal')
	for e in ells:
	    ax.add_artist(e)
	    e.set_clip_box(ax.bbox)
	    e.set_facecolor('blue')
	#ax.imshow(slices['ax1'], origin='lower', cmap='hsv', interpolation='bilinear')
	plt.show()


if __name__ == "__main__":
	optparser = optparse.OptionParser()
	optparser.add_option(
	    "-s", "--sensor_data", default="intel.script",
	    help="Sensor data file location"
	)
	opts = optparser.parse_args()[0]
	maplocation = opts.sensor_data
	with open(maplocation, 'r') as file:
		log = file.read().split('\n')

	csm = CSM()
	state = [0.0, 0.0, 0.0]
	prev_scan = 0.0
	for line in log[425:432]:
		if line.startswith('POS'):
			prev_state = state
			state = [float(val) for val in line.split()[3:]]
		elif line.startswith('LASER-RANGE'):
			tok = line.split()
			minRange = convert_to_radians(float(tok[3])-90)
			maxRange = convert_to_radians(float(tok[5][:-1])-90)
			numScans = int(tok[4])
			increment = (maxRange-minRange)/numScans
			if prev_scan != 0.0:
				#print "Building lookup table"
				csm.buildLookupTable(prev_state, prev_scan)
			else:
				scan = tok[6:]
				prev_scan = scan
				continue
			scan = tok[6:]
			prev_scan = scan
			rad = convert_to_radians(15)
			slice = np.zeros(shape=(100,100))
			i = 0
			#slices = {}
			slices = np.zeros(shape=(30, 100, 100))
			orientation = state[2]
			for theta in np.arange(-rad+state[2], rad+state[2], 2*rad/30):
				state[2] = theta
				#print state
				prev_slice = slice
				slice = csm.twoDslices(scan, state)
				slice = slice*100
				slices[i] = slice
				i+=1
				if (prev_slice == slice).all():
					print "THEY ARE THE SAME FOR SLICE CALCULATION!!"
			index = slices.argmax()
			index_tuple = np.unravel_index(index, slices.shape)
			print index_tuple
			slice_num = index_tuple[0]
			plotGaussian(slice_num, state)