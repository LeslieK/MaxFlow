'''
Module includes classes to process flow networks.

Classes:
FlowEdge, FlowNetwork

'''
from collections import Counter
import random

class FlowEdge(object):
	'''an edge is a flow network, flow/capacity'''

	def __init__(self, v, w, capacity, flow=None):
		'''create a flow edge v -> w'''
		self._v = v
		self._w = w
		self._capacity = capacity
		if flow is not None:
			self._flow = flow
		else:
			self._flow = 0

	def source(self):
		'''the vertex from which the edge leaves'''
		return self._v

	def sink(self):
		'''the vertex on which the edge terminates'''
		return self._w

	def other(self, vertex):
		'''the vertex at the other end of the edge'''
		if vertex == self._v:
			return self._w
		else:
			return self._v

	def capacity(self):
		'''capacity of this edge'''
		return self._capacity

	def flow(self):
		'''flow of this edge'''
		return self._flow

	def residualCapacityTo(self, vertex):
		'''the residual capacity toward vertex'''
		if vertex == self._v:
			# backward edge
			return self._flow
		else:
			# forward edge
			return self._capacity - self._flow

	def addResidualFlowTo(self, vertex, delta):
		'''add delta flow toward v'''
		if vertex == self._v:
			# backward edge
			self._flow -= delta
		else:
			# forward edge
			self._flow += delta

	def __repr__(self):
		'''unambiguous description of flow edge'''
		return '{} -> {}  {}/{}'.format(self._v, self._w, self._flow, self._capacity)

class FlowNetwork(object):
	'''graph representing flow network'''
	def __init__(self, filename=None, V=None, E=None):
		'''create an empty flow network with V vertices'''
		if filename is not None:
			# read from file
			with open(filename, 'r') as f:
				V = int(f.readline())
				E = int(f.readline())
				self._V = V
				self._E = 0
				self._adj = []
				for v in range(V):
					self._adj.append(Counter())
				for line in f.readlines():
					terms = line.split()
					v = int(terms[0])
					w = int(terms[1])
					capacity = int(float((terms[2])))
					self.addEdge(FlowEdge(v, w, capacity))
				return
		# no filename
		self._V = V
		self._E = 0
		self._adj = []
		for v in range(V):
			self._adj.append(Counter())
		if E is not None:
			# generate a random network
			for _ in range(E):
				v = int(round(random.uniform(0, V-1)))
				w = int(round(random.uniform(0, V-1)))
				capacity = int(round(random.uniform(1, 100)))
				self.addEdge(FlowEdge(v, w, capacity))
				self._E += 1

	def addEdge(self, edge):
		'''add flow edge to this flow network'''
		v = edge.source()
		w = edge.sink()
		self._adj[v][edge] += 1
		self._adj[w][edge] += 1
		self._E += 1

	def adj(self, v):
		'''an iterable of all edges incident on v'''
		return self._adj[v]

	def edges(self):
		'''an iterable of all edges in flow network'''
		edges = []
		for v in range(self._V):
			for e in self.adj(v):
				if e.sink() != v:
					# do not include self-loops
					edges.append(e)
		return edges

	def V(self):
		'''return the number of vertices'''
		return self._V

	def E(self):
		'''return the number of edges'''
		return self._E

	# def __repr__:
	# 	'''unambiguously describes the flow network'''
	# 	pass
		

