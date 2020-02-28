class UnionFind:
	"""Weighted quick-union with path compression.
	The original Java implementation is introduced at
	https://www.cs.princeton.edu/~rs/AlgsDS07/01UnionFind.pdf
	"""

	def __init__(self, n):
		self._id = list(range(n))
		self._sz = [1] * n
	def find(self, i):
		j = i
		while (j != self._id[j]):
			self._id[j] = self._id[self._id[j]]
			j = self._id[j]
		return j
	def union(self, p, q):
		i = self.find(p)
		j = self.find(q)
		if i == j: return

		if (self._sz[i] < self._sz[j]):
			self._id[i] = j
			self._sz[j] += self._sz[i]
		else:
			self._id[j] = i
			self._sz[i] += self._sz[j]
	def __getitem__(self,i):
		return self.find(i)