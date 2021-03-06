import collections
from decimal import Decimal
_INF = Decimal('infinity')


class FF(object):
	"""Ford-Fulkerson algorithm for computing maxflow/mincut using shortest augmenting path rule."""

	def __init__(self, G, s, t):
		'''computes maxflow/mincut in flow network G from source s to target t'''

		if s < 0 or s > G.V():
			raise ValueError('s is out of range: {}'.format(s))
		if t < 0 or t > G.V():
			raise ValueError('t is out of range: {}'.format(t))
		if s == t:
			raise ValueError('bad input: s equals t')
		# current value of flow (excess at target)
		self._value = self._excess(G, t)

		if (not self._isFeasible(G, s, t)):
			raise ValueError('initial flow is infeasible; violates network constraints')
		# initialize data structures (size of all connected components)
		self._marked = [False] * G.V()
		self._edgeTo = [-1] * G.V()

		# while there exists an augmenting path, add flow to it
		while(self._hasAugmentingPath(G, s, t)):
			# find bottleneck
			bottle = _INF
			v = t
			while v != s:
				e = self._edgeTo[v]
				bottle = min(e.residualCapacityTo(v), bottle)
				v = e.other(v)
			# increment flow on augmenting path by amount bottle
			v = t
			while v != s:
				e = self._edgeTo[v]
				e.addResidualFlowTo(v, bottle)
				v = e.other(v)
			# update flow
			self._value += bottle

	def _hasAugmentingPath(self, G, s, t):
		'''is there an augmenting path? if yes, store path in edgeTo'''

		# this algorithm uses shortest path (BFS, queue) to search for augmenting paths
		# alternatives: fattest path (priority queue), DFS path ('stack')

		# BFS (breadth-first-search)
		self._marked = [False] * G.V()
		self._edgeTo = [-1] * G.V()

		q = collections.deque()
		q.append(s)
		self._marked[s] = True
		while len(q) > 0:
			# q is not empty
			v = q.popleft()
			for e in G.adj(v):
				w = e.other(v)
				if e.residualCapacityTo(w) > 0:
					if not self._marked[w]:
						self._marked[w] = True
						self._edgeTo[w] = e
						q.append(w)
		return self._marked[t]

	def value(self):
		'''return value of max flow'''
		return self._value

	def inCut(self, v):
		'''is v in the s side of the min s-t cut?'''
		return self._marked[v]

	def _excess(self, G, v):
		'''return net flow at vertex v'''
		excess = 0.0
		for e in G.adj(v):
			if v == e.source():
				# edge leaves v
				excess -= e.flow()
			else:
				# edge points to v
				excess += e.flow()
		return excess

	def _isFeasible(self, G, s, t):
		EPSILON = 10e-11
		# check flow constraints are met
		#for v in range(G.V()):
		for v in range(G.V()):  # vertices 
			for e in G.adj(v):
				if e.flow() < 0: 
					print 'e.flow < 0'
					return False
				elif e.flow() > e.capacity():
					print 'e.flow() > e.capacity()'
					return False
				if e.flow() < 0 or (e.flow() > e.capacity()):
				 	print 'Edge does not satisfy constraints'
				 	return False

		# check source and target
		if abs(self._value + self._excess(G, s)) > EPSILON:
			print "excess flow at source"
			return False
		if abs(self._value - self._excess(G, t)) > EPSILON:
			print "excess flow at target"
			return False
		# check that net flow into each vertex is 0, except for s and t
		#for v in range(G.V()):
		for v in range(G.V()):
			if v == s or v == t:
				continue
			elif abs(self._excess(G, v)) > EPSILON:
				print 'Net flow out of v does not equal zero'
				return False
		return True






