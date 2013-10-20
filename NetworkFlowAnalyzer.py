'''
Performs network flow analysis on a network of V points

Each point is a bike station in the citibike network
Pairs of points are connected if they are a distance of _DISTANCE km
from each other.
Capacity on an edge is the sum of the bike docks at each edge vertex.

usage: run NetworkFlowAnalyzer "citybike.json" <start-station> <end-station>

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
parser.add_argument("--stations", action="store_true", help="a list of stations from input file")
parser.add_argument("--start_station", "-s", default="Broadway & Battery Pl", help="a station name from input file")
parser.add_argument("--end_station", "-e", default="E 47 St & 1 Ave", help="a station name from input file", type=str)
args = parser.parse_args()

# read network data
with open(args.filename) as f:
    input = json.load(f)

num_stations = len(input["stationBeanList"])
stations = np.random.random((num_stations, 2))  # 2 cols for latitude, longitude
#totalDocks = np.zeros((num_stations,), dtype=np.uint8)
totalDocks = np.zeros((num_stations,), dtype=np.int16)
names = {}  # maps vertex number to stationNam
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

# process optional args
if args.stations:
    # print list of stations and then return to command line
    for row in range(num_stations):
        print names[row]
else:
    # strategy:
    start_station = args.start_station
    end_station = args.end_station
    # get vertex numbers
    source = vertex[start_station]
    target = vertex[end_station]
    # connect stations
    flow_edges = MaxFlowUtils.findFlowEdges(stations, totalDocks, source, MaxFlowUtils.getCapacities_Max)

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
        # plot flow network
        MaxFlowUtils.plotFlowNetwork(stations, source, target, flow_edges)
        plt.show()
        # flow_edges.append(FlowEdge(329, 330, totalDocks[329] + totalDocks[330]))

        # plot flow
        flowpath = MaxFlowUtils.flowPath(stations, flownet, source)
        flowpath2 = copy.copy(flowpath)
        MaxFlowUtils.plotFlow(stations, source, target, flowpath)
        plt.show()

        # convert vertices in flow path to station names
        station_path = MaxFlowUtils.toStationNames(flowpath2, names)
        while len(station_path) > 0:
            print station_path.popleft()

        print 'start: {}'.format(start_station)
        print 'end: {}'.format(end_station)
        print 'maxflow = {}'.format(maxflow.value())
        print

        mincut = MaxFlowUtils.findMinCut(maxflow, num_stations)
        stcut, stcut_v = MaxFlowUtils.findSTcut(mincut, flownet, names)
        print 'Edges that, if cut, would separate {} from {} (aka st-cut):\n'.format(start_station, end_station)
        for e in stcut:
            print e

        # plot edges in st-cut
        MaxFlowUtils.plotSTcut(stations, source, target, flowpath, stcut_v)
        plt.show()
    else:
        print '{} and {} are not connected.'.format(start_station, end_station)

    




