'''
Performs network flow analysis on a network of V points

Each point is a bike station in the citibike network
Pairs of points are connected if they are a distance of _DISTANCE km
from each other.
Capacity on an edge is calculated from the bike docks at each edge vertex.

usage: run NetworkFlowAnalyzer "citybike.json" <start-station> <end-station>

'''
import json
import argparse
import matplotlib.pyplot as plt
import numpy as np
from MaxFlowUtils import getCapacities_Max  
from FordFulkerson import FF
#from FattestPath import FF
from Citibike import Station, Graph

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
stations = np.zeros((num_stations, 3))
totalDocks = np.zeros((num_stations,), dtype=np.int16)  # need to support neg numbers
nodesByName = {}        # maps stationName => node
nodesByNumber = {}     # maps node number => node object

for row, station in enumerate(input["stationBeanList"]):
    stations[row][0] = station["latitude"]
    stations[row][1] = station["longitude"]
    stations[row][2] = station["totalDocks"]
    node = Station(station, row)
    nodesByNumber[node.number] = node
    if node.name == args.start_station:
        start_station = node
    elif node.name == args.end_station:
        end_station = node


# process optional args
if args.stations:
    # print list of stations and then return to command line
    for row in range(num_stations):
        print nodesByNumber[row].name
else:
    graph = Graph(stations, start_station, end_station, getCapacities_Max, nodesByNumber)
    
    if graph.isConnected:
        # plot flow network
        graph.plotFlowNetwork()

        # plot flow
        graph.plotFlow()
        plt.show()

        # plot S-T cut
        graph.plotSTcut()
        plt.show()

        #print flow
        graph.printFlow(nodesByNumber)

        print 'start: {}'.format(start_station.name)
        print 'end: {}'.format(end_station.name)
        print 'maxflow = {}'.format(graph.maxflow())
        print

        graph.printSTcut(nodesByNumber)
    else:
        print '{} and {} are not connected.'.format(start_station.name, end_station.name)