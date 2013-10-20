import math, collections
from MaxFlowLib import FlowEdge
import matplotlib.pyplot as plt
import numpy as np
import copy
_DISTANCE = .4 	# km

# vectorized implementation
def distance(v, stations):
	'''find distance from v to all other vertices'''
	R = 6371		
	dLat_dLong = np.radians(stations[:, [0, 1]] - stations[v, [0, 1]]) 
	dLat = dLat_dLong[:, 0]
	dLong = dLat_dLong[:, 1]
	lat = np.radians(stations[:, 0])  		# col vector of lat values for all stations; num_stations x 1
	lat_v = np.radians(stations[v][0])    	# scalar for this station

	dLatH = np.divide(dLat, 2)
	dLongH =  np.divide(dLong, 2)
	a = np.sin(dLatH) ** 2 + np.sin(dLongH) ** 2 * np.cos(lat) * np.cos(lat_v)
	c = np.arctan2(np.sqrt(a), np.sqrt(1-a)) * 2
	d = c * R
	return d

# vectorized implementation that calls distance2
def findFlowEdges(stations, totalDocks, source, capacity_func):
	'''return flow edges

	stations = numpy array; column 0: latitude  column 1: longitude  
	totalDocks = numpy array: column 0 is totalDocks at station 
	source = source vertex
	'''
	flow_edges = collections.deque()
	num_stations = len(stations)
	marked_v = np.ones((num_stations,), dtype=bool) # marks all from vertices
	marked_w = np.ones((num_stations,), dtype=bool) # marks all to vertices
	q = collections.deque()
	q.append(source)
	
	while len(q) > 0:
		# q not empty
		v = q.popleft()
		marked_v[v] = False
		d = distance(v, stations)
		connected = np.nonzero(d < _DISTANCE)[0]
		capacities = capacity_func(totalDocks[connected], totalDocks[v])
		# filter connected vertices; exclude if marked_v[w] = 0
		filter_v = marked_v[connected]
		capacities = capacities[filter_v]
		connected = connected[filter_v]
		# only 2 for loops!!
		for i, w in enumerate(connected):
			flow_edges.append(FlowEdge(v, w, capacities[i]))
		# filter w before appending to q
		filter_w = marked_w[connected]
		connected = connected[filter_w]
		for w in connected:
			q.append(w)
			marked_w[w] = False
	return flow_edges

def getCapacities_Max(array_of_capacities, capacity_of_v):
	'''capacity of edge = max of totalDocks at v, w

	return array of capacities'''
	capacities = array_of_capacities
	delta = capacities - capacity_of_v
	capacities[delta < 0] = capacity_of_v
	return capacities

def getCapacities_Sum(array_of_capacities, capacity_of_v):
	'''capacity of edge = sum of totalDocks at v, w

	return array of capacities'''
	return array_of_capacities + capacity_of_v

def drawScatterPlot(stations, source, target):
	'''draw scatter plot

	stations: num_stations x 2 numpy array
	source: source vertex
	target: target vertex
	'''
	# draw scatter plot
	plt.xlabel('latitude')
	plt.ylabel('longitude')
	plt.title('citibike stations in NYC')
	plt.plot(stations[:, 0], stations[:, 1], 'ro')

	# annotate start, target stations
	plt.annotate('start', xy=(stations[source][0], stations[source][1]), xycoords='data', xytext=(-50,+50), \
	textcoords='offset points', arrowprops=dict(arrowstyle="->", connectionstyle="arc3", linewidth=3.5), fontsize=16)
	plt.annotate('end', xy=(stations[target][0], stations[target][1]), xycoords='data', xytext=(-50, +50), \
	textcoords='offset points',  arrowprops=dict(arrowstyle="->", connectionstyle="arc3", linewidth=3.5), fontsize=16)


def plotFlow(stations, source, target, flowpath):
	'''plot path of flow'''

	drawScatterPlot(stations, source, target)
	plt.title('citibike stations in NYC: path of flow')

	# draw path of flow
	flowpath2 = copy.copy(flowpath)
	while len(flowpath2) > 0:
		e = flowpath2.popleft()
		v = e.source()
		w = e.sink()
		plt.plot([stations[v][0],stations[w][0]], [stations[v][1], stations[w][1]], 'b-')
	#plt.show()


def plotFlowNetwork(stations, source, target, flow_edges):
	'''draw network of NYC citibike stations

	stations: num_stations x 2 numpy array
	flow_edges: a list
	source: source vertex
	target: target vertex
	'''
	drawScatterPlot(stations, source, target)
	plt.title('citibike stations in NYC: flow network')

	# draw lines
	for e in range(len(flow_edges)):
	    f = flow_edges[e]
	    v = f.source()
	    w = f.sink()
	    plt.plot([stations[v][0],stations[w][0]], [stations[v][1], stations[w][1]], 'b-')
	#plt.show()

def plotSTcut(stations, source, target, flowpath, stcut_v):
	'''plot mincut and the edges that bridge the st-cut

	stcut_v: set of edges that bridge st-cut
	'''
	plotFlow(stations, source, target, flowpath)
	plt.title('citibike stations in NYC: edges in st-cut')

	# draw edges in st-cut
	for e in stcut_v:
		v = e.source()
		w = e.sink()
		plt.plot([stations[v][0], stations[w][0]], [stations[v][1], stations[w][1]], color="black", linewidth=4.0)
	#plt.show()


def findMinCut(maxflow, num_stations):
	'''
	find all stations in mincut 
	'''
	mincut = set()
	for v in range(num_stations):
		if maxflow.inCut(v):
			mincut.add(v)
	return mincut

def findSTcut(mincut, flownet, names):
	'''prints all edges in s-t cut

	s-t cut: set of edges from partition containing s
	to partition containing t
	'''
	stcut = set()		# set of edges using station names
	stcut_v = set()  	# set of edges using vertices
	#tot = 0
	for vertex in mincut:
		for e in flownet.adj(vertex):
			if vertex == e.source() and e.other(vertex) not in mincut: 
				if e.flow() == e.capacity():
					#tot += e.capacity()
					stcut.add('{:<30} {:^5} {:<30}  {:>5}/{:<5}'.format(names[vertex], '=>', names[e.sink()], e.flow(), e.capacity()))
					stcut_v.add(e)
	#print tot
	return (stcut, stcut_v)

def flowPath(stations, flownet, source):
	'''return vertices in flow path

	flownet: flow network
	source = source vertex
	'''
	# build flow path
	flowpath = collections.deque()
	# working queue to find next edge in path
	q = collections.deque()
	num_stations = len(stations)
	marked = [False] * num_stations
	marked[source] = True
	q.append(source)
	while len(q) > 0:
		v = q.popleft()
		for e in flownet.adj(v):
			if v == e.source() and e.flow() > 0:
				flowpath.append(e)
				if not marked[e.sink()]:
					marked[e.sink()] = True
					q.append(e.sink())
	return flowpath

def toStationNames(flowpath, names):
	'''convert each flow edge to station names

	flowpath: deque of flow edges
	'''
	station_path = collections.deque()
	while len(flowpath) > 0:
		e = flowpath.popleft()
		source_station = names[e.source()]
		target_station = names[e.sink()]
		edge = '{:<30} {:^5} {:<30}  {:>5}/{:<5}\n'.format(source_station, '=>', target_station, e.flow(), e.capacity())
		station_path.append(edge)
	return station_path






