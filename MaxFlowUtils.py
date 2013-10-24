import math, collections
from MaxFlowLib import FlowEdge
import matplotlib.pyplot as plt
import numpy as np
import copy
_DISTANCE = .4 	# km

def getCapacities_Max(array_of_capacities, this_station_capacity):
	'''capacity of edge = max of totalDocks at v, w

	v, w: stations at each end of edge
	v = from_station
	w = to_station
	return array of capacities'''
	# return [max(connected_station, this_station_capacity) for connected_station in array_of_capacities]
	delta = array_of_capacities - this_station_capacity
	array_of_capacities[delta < 0] = this_station_capacity
	return array_of_capacities

def getCapacities_Sum(array_of_capacities, this_station_capacity):
	'''capacity of edge = sum of totalDocks at v, w

	v, w: stations at each end of edge
	v = from_station
	w = to_station
	return array of capacities'''
	return array_of_capacities + this_station_capacity

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









