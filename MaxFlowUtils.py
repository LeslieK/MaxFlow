import math, collections
from MaxFlowLib import FlowEdge
import matplotlib.pyplot as plt
import numpy as np
import copy
_DISTANCE = .4 	# km

def distance(lat1, long1, lat2, long2):
	'''returns distance, in km, between 2 points
	'''
	R = 6371 # Earth's radius in km
	dLat = math.radians(lat2-lat1)
	dLon = math.radians(long2-long1)
	lat1 = math.radians(lat1)
	lat2 = math.radians(lat2)

	a = (math.sin(dLat/2) * math.sin(dLat/2) +
	    math.sin(dLon/2) * math.sin(dLon/2) * math.cos(lat1) * math.cos(lat2))
	c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a)) 
	d = R * c
	return d

# for row in range(num_stations):
# 	# dLat = math.radians(lat2-lat1) => dLat_dLong[:, 0]		num_stations x 1
# 	# dLon = math.radians(long2-long1) => dLat_dLong[:, 1]		num_stations x 1
# 	dLat_dLong = stations[:, [0, 1]] - stations[row, [0, 1]] #  num_stations x 2
# 	dLat = map(math.radians, dLat_dLong[:, 0])
# 	dLon = map(math.radians, dLat_dLong[:, 1])
# 	lat1 = map(math.radians, stations[:, 0])  # col vector of lat1 values for all stations; num_stations x 1
# 	lat2 = math.radians(stations[row][0])     # scalar for this station (this row)

# 	a = (math.sin(dLat/2) * math.sin(dLat/2) +
# 	    math.sin(dLon/2) * math.sin(dLon/2) * math.cos(lat1) * math.cos(lat2))

# dLatH = map(lambda x: x/2., dLat)
# dLongH =  map(lambda x: x/2., dLong)
# np.sin(dLatH) ** 2 + np.sin(dLonH) ** 2 * np.cos(lat1) * np.cos(lat2)

def findFlowEdges(stations, totalDocks, source):
	'''return flow edges

	stations = numpy array; column 0: latitude  column 1: longitude  
	capacities[row] = capacity of vertex=row  
	source = source vertex
	'''
	flow_edges = collections.deque()

	num_stations = len(stations)
	marked_v = [False] * num_stations # marks all from vertices
	marked_w = [False] * num_stations # marks all to vertices
	q = collections.deque()
	q.append(source)
	
	while len(q) > 0:
		# q not empty
		v = q.popleft()
		marked_v[v] = True
		vlat = stations[v][0]
		vlong = stations[v][1]
		for w in range(num_stations):
			if not marked_v[w]:
				wlat = stations[w][0]
				wlong = stations[w][1]
				d = distance(vlat, vlong, wlat, wlong)
				if d < _DISTANCE:
					capacity = np.max(totalDocks[v] + totalDocks[w])
					# first edge is from source; last edge is incident on target
					flow_edges.append(FlowEdge(v, w, capacity)) 
					if not marked_w[w]:
						# w has not been put on queue
						marked_w[w] = True
						q.append(w)
	return flow_edges

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






