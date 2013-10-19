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
from MaxFlowUtils import _DISTANCE  # km
from MaxFlowLib import FlowEdge, FlowNetwork
from FordFulkerson import FF
import ConnectedComponent
import copy
from collections import defaultdict


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
start_station = "Broadway & Battery Pl"
end_station = "E 47 St & 1 Ave"
# get vertex numbers
source = vertex[start_station]
target = vertex[end_station]
# connect stations
flow_edges = MaxFlowUtils.findFlowEdges(stations, totalDocks, source)

# build flow network
# number of vertices = num_stations
flownet = FlowNetwork(V=num_stations)
for e in flow_edges:
	flownet.addEdge(e)

# get number of connected components in graph
cc = ConnectedComponent.CC(flownet)
if cc.id(source) == cc.id(target):
    # calculate maxflow using Ford-Fulkerson algorithm
    maxflow = FF(flownet, source, target)
    # draw scatter plot
    MaxFlowUtils.drawScatterPlot(stations, flow_edges, source, target)
    # flow_edges.append(FlowEdge(329, 330, totalDocks[329] + totalDocks[330]))

    # plot flow path
    flowpath = MaxFlowUtils.flowPath(stations, flownet, source)
    flowpath2 = copy.copy(flowpath)
    MaxFlowUtils.plotFlow(stations, source, target, flowpath)
    # convert vertices in flow path to station names
    station_path = MaxFlowUtils.toStationNames(flowpath2, names)
    while len(station_path) > 0:
        print station_path.popleft()

    print 'start: {}'.format(start_station)
    print 'target: {}'.format(end_station)
    print 'maxflow = {}'.format(maxflow.value())
else:
    print '{} and {} are not connected.'.format(start_station, end_station)

# mincut = MaxFlowUtils.findMinCut(maxflow, num_stations)
# for v in mincut:
#     print names[v] + "\n"




