'''
Performs network flow analysis on a network of V points
Each point is a bike station in the citibike network
Pairs of points are connected if they are a distance of _DISTANCE km
from each other.
Capacity on an edge is the sum of the bike docks at each edge vertex.
'''
import json
import argparse
import matplotlib.pyplot as plt
import numpy as np
import random
import MaxFlowUtils
from MaxFlowLib import FlowEdge, FlowNetwork
from FordFulkerson import FF
import ConnectedComponent

from collections import defaultdict
_DISTANCE = .4  # km

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="filename of json network", type=str)
args = parser.parse_args()

# read network data
with open(args.filename) as f:
    input = json.load(f)

num_stations = len(input["stationBeanList"])
stations = np.random.random((num_stations, 2))  # 2 cols for latitude, longitude
totalDocks = np.zeros((num_stations,), dtype=np.uint8)
names = {}  # maps vertex number to stationName
vertex = {} # maps stationName to vertex number

# store data
row = 0
for station in input["stationBeanList"]:
	# store data in numpy arrays
    totalDocks[row] = station["totalDocks"]
    stations[row][0] = station["latitude"]
    stations[row][1] = station["longitude"]
    vertex[station["stationName"]] = row
    names[row] = station["stationName"]
    row += 1

# strategy:
start_station = "1 Ave & E 15 St"
end_station = "Broadway & W 51 St"
# build a digraph from start to end
source = vertex[start_station]
target = vertex[end_station]
# connect stations
flow_edges = MaxFlowUtils.findFlowEdges(stations, totalDocks, source, target)

# flow_edges = []
# for v in range(num_stations):
# 	vlat = stations[v][0]
# 	vlong = stations[v][1]
# 	for w in range(num_stations):
# 		if w > v:
# 			wlat = stations[w][0]
# 			wlong = stations[w][1]
# 			d = MaxFlowUtils.distance(vlat, vlong, wlat, wlong)
# 			if d < _DISTANCE:
# 				capacity = totalDocks[v] + totalDocks[w]
# 				flow_edges.append(FlowEdge(v, w, capacity))
#                 #flow_edges.append(FlowEdge(w, v, capacity))
# flow_edges.append(FlowEdge(329, 330, totalDocks[329] + totalDocks[330]))

#draw scatter plot
plt.xlabel('latitude')
plt.ylabel('longitude')
plt.title('citibike stations in NYC')
plt.plot(stations[:, 0], stations[:, 1], 'ro')
# draw lines
lat_array = np.arange(len(flow_edges) * 2)
long_array = np.arange(len(flow_edges) * 2)
for e in range(len(flow_edges)):
    f = flow_edges[e]
    v = f.source()
    w = f.sink()
    plt.plot([stations[v][0],stations[w][0]], [stations[v][1], stations[w][1]], 'b-')
plt.show()

# build flow network
# number of vertices = num_stations
flownet = FlowNetwork(num_stations)
for e in flow_edges:
	flownet.addEdge(e)

# # find connected components in flow network
# cc = ConnectedComponent.CC(flownet)
# print 'Number of connected components: {}'.format(cc.count())
# counts = [0] * cc.count()
# vertices = defaultdict(list)
# for v in range(num_stations):
#     counts[cc.id(v)] += 1
#     vertices[cc.id(v)].append(v)

# # run Ford-Fulkerson algorithm over graph to find max flow
# start = vertex["1 Ave & E 15 St"]
# end = vertex["Broadway & W 51 St"]

# # check if connected
# if cc.id(start) == cc.id(end):
#     list_of_vertices =  vertices[cc.id(start)]
#     maxflow = FF(flownet, start, end, list_of_vertices)
# else:
# 	print '{} is not connected to {}'.format(names[start], names[end])



