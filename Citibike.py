from MaxFlowLib import FlowEdge, FlowNetwork
from FordFulkerson import FF
#from FattestPath import FF
import ConnectedComponent
import matplotlib.pyplot as plt
import numpy as np
import collections
from MaxFlowUtils import getCapacities_Max, drawScatterPlot
_DISTANCE = .4 	# km

class Station(object):
	'''
	encapsulates a citibike station
	'''
	def __init__(self, station, node_number):
		# store data in numpy arrays
		self.number = node_number
		self.latitude = station["latitude"]
		self.longitude = station["longitude"]
		self.totalDocks = station["totalDocks"]
		self.name = station["stationName"]

	# vectorized implementation
	def distance(self, stations):
		'''find distance from this station to all other stations'''
		R = 6371		
		dLat_dLong = np.radians(stations[:, [0, 1]] - stations[self.number, [0, 1]]) 
		dLat = dLat_dLong[:, 0]
		dLong = dLat_dLong[:, 1]
		lat = np.radians(stations[:, 0])  		# col vector of lat values for all stations; num_stations x 1
		lat_v = np.radians(stations[self.number][0])    	# scalar for this station

		dLatH = np.divide(dLat, 2)
		dLongH =  np.divide(dLong, 2)
		a = np.sin(dLatH) ** 2 + np.sin(dLongH) ** 2 * np.cos(lat) * np.cos(lat_v)
		c = np.arctan2(np.sqrt(a), np.sqrt(1-a)) * 2
		d = c * R
		return d

class Graph(object):
	'''build network of stations; connect with edges'''

	def __init__(self, stations, start_station, end_station, capacity_func, nodesByNumber):
		'''build graph of stations and edges

		start_station: station object
		end_station: station_object
		'''

		# global
		# self.nodesByName = nodesByName
		# self.nodesByNumber = nodesByNumber
		
		# flow network
		self.stations = stations
		self.num_stations = len(stations)
		self.flownet = FlowNetwork(V=self.num_stations)
		self.flowpath = []

		# start and end nodes
		self.start_station = start_station
		self.end_station = end_station
		self.isConnected = False

		# edges
		flow_edges = set()

		# find edges in flow network
		marked_v = np.ones((self.num_stations,), dtype=bool) # marks all from_vertices (aka v)
		marked_w = np.ones((self.num_stations,), dtype=bool) # marks all to_vertices (aka w)
		
		q = collections.deque()			# a queue of station objects
		q.append(self.start_station)
	
		while len(q) > 0:
			# q not empty
			station = q.popleft()

			# from station v
			v = station.number
			marked_v[v] = False
			# d is an array of distances
			d = station.distance(self.stations)
			# connected is array of station numbers that connect to this station
			connected = np.nonzero(d < _DISTANCE)[0]
			capacities = capacity_func(self.stations[:, 2][connected], station.totalDocks)
			# filter connected vertices; exclude if marked_v[w] = 0
			filter_v = marked_v[connected]
			capacities = capacities[filter_v]
			connected = connected[filter_v]
			# only 2 for loops!!
			for i, w in enumerate(connected):
				flow_edges.add(FlowEdge(v, w, capacities[i]))
			# filter w before appending to q; exclude all w that are marked
			filter_w = marked_w[connected]
			connected = connected[filter_w]
			for w in connected:
				# w is node number
				node_w = nodesByNumber[w]
				q.append(node_w)
				marked_w[w] = False

		# add edges to network
		for e in flow_edges:
			self.flownet.addEdge(e)
		# check if start and end are connected
		self.setIsConnected()
		if self.isConnected:
			self.ff = FF(self.flownet, self.start_station.number, self.end_station.number)
			self.flowpath = self.flowPath()

	def setIsConnected(self):
		'''sets the instance variable isConnected

		True if graph is one connected component; False otherwise
		'''
		cc = ConnectedComponent.CC(self.flownet)
		source = self.start_station.number
		end = self.end_station.number
		self.isConnected = (cc.id(source) == cc.id(end))

	def maxflow(self):
		'''find maxflow of flow network'''
		return self.ff.value()


	def plotFlowNetwork(self):
		'''draw network of NYC citibike stations

		source: start_station number
		end: end_station number
		'''
		source = self.start_station.number
		end = self.end_station.number
		drawScatterPlot(self.stations, source, end)
		plt.title('citibike stations in NYC: flow network')

		# draw lines
		for f in self.flownet.edges():
		    v = f.source()
		    w = f.sink()
		    plt.plot([self.stations[v][0],self.stations[w][0]], [self.stations[v][1], self.stations[w][1]], 'b-')
		plt.show()

	def plotFlow(self):
		'''plot path of flow'''
		source = self.start_station.number
		end = self.end_station.number
		drawScatterPlot(self.stations, source, end)
		plt.title('citibike stations in NYC: paths of flow')

		for e in self.flowpath:
			v = e.source()
			w = e.sink()
			plt.plot([self.stations[v][0], self.stations[w][0]], [self.stations[v][1], self.stations[w][1]], 'b-')


	def printFlow(self, nodesByNumber):
		'''print each flow edge using station names'''
		station_path = self.toStationNames(self.flowpath, nodesByNumber)
		while len(station_path) > 0:
			print station_path.popleft()


	def flowPath(self):
		'''return a list of edges in flow path

		flownet: flow network
		source = start_station number
		'''
		source = self.start_station.number
		# build flow path
		flowpath = []
		# working queue to find next edge in path
		q = collections.deque()
		marked = [False] * self.num_stations
		marked[source] = True
		q.append(source)
		while len(q) > 0:
			v = q.popleft()
			for e in self.flownet.adj(v):
				if v == e.source() and e.flow() > 0:
					flowpath.append(e)
					if not marked[e.sink()]:
						marked[e.sink()] = True
						q.append(e.sink())
		return flowpath


	def toStationNames(self, flowpath, nodesByNumber):
		'''convert each flow edge to station names

		flowpath: list of edges
		'''
		station_path = collections.deque()
		for e in self.flowpath:
			# source_station = self.nodesByNumber[e.source()].name
			# end_station = self.nodesByNumber[e.sink()].name
			source_station = nodesByNumber[e.source()].name
			end_station = nodesByNumber[e.sink()].name
			edge = '{:<30} {:^5} {:<30}  {:>5}/{:<5}\n'.format(source_station, '=>', end_station, e.flow(), e.capacity())
			station_path.append(edge)
		return station_path

	def findMinCut(self):
		'''
		find all stations in mincut 
		'''
		mincut = set()
		for v in range(self.num_stations):
			if self.ff.inCut(v):
				mincut.add(v)
		return mincut

	def findSTcut(self, nodesByNumber=None, names=False):
		'''
		find all edges that cross the s-t cut

		s-t cut: a partition of the flow network
		mincut: all stations connected to the start_station after no more flow edges can be found
		'''
		mincut = self.findMinCut()
		stcut = set()		# set of edges using station names
		for vertex in mincut:
			for e in self.flownet.adj(vertex):
				if vertex == e.source() and e.other(vertex) not in mincut: 
					if e.flow() == e.capacity():
						if names:
							# from_station = self.nodesByNumber[vertex].name
							# to_station = self.nodesByNumber[e.sink()].name
							from_station = nodesByNumber[vertex].name
							to_station = nodesByNumber[e.sink()].name
							stcut.add('{:<30} {:^5} {:<30}  {:>5}/{:<5}'.format(from_station, '=>', to_station, e.flow(), e.capacity()))
						else:
							stcut.add(e)
		return stcut

	def printSTcut(self, nodesByNumber):
		'''print all edges that cross the s-t cut

		s-t cut: a partition of the flow network
		'''
		stcut = self.findSTcut(nodesByNumber, names=True)
		print 'Edges that, if cut, would separate {} from {} (aka st-cut):\n'.format(self.start_station.name, self.end_station.name)
		for e in stcut:
			print e


	def plotSTcut(self):
		'''plot mincut and the edges that bridge the st-cut

		stcut: set of edges that bridge st-cut
		'''
		stcut = self.findSTcut(names=False)
		self.plotFlow()
		plt.title('citibike stations in NYC: edges in st-cut')

		# draw edges in st-cut
		for e in stcut:
			v = e.source()
			w = e.sink()
			plt.plot([self.stations[v][0], self.stations[w][0]], [self.stations[v][1], self.stations[w][1]], color="black", linewidth=4.0)
		#plt.show()


