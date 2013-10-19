'''
Performs network flow analysis on a network of V points
Each point is a bike station in the citibike network
Pairs of points are connected if they are a distance of _DISTANCE km
from each other.
Capacity on an edge is the sum of the bike docks at each edge vertex.
'''

from MaxFlowLib import FlowNetwork
from FordFulkerson import FF

# Flownet(filename) or Flownet(V=n) or Flownet(V=n, E=m)
flownet = FlowNetwork(V=100, E=300)
V = flownet.V()

maxflow = FF(flownet, 0, V-1)

print maxflow.value()
flownet.edges()

