import numpy as np
from numpy.random import random


class Network:
    """
    Abstract a network graph.
    """

    @staticmethod
    def AlbertBarabasi(n, directed=False, initial_graph=None):
        if initial_graph is None:
            # The bootstrap graph should have a minimum of 2 connected nodes
            g = Network.FullyConnected(2, directed, allow_diagonals=False)
            first_node_ix = 2
        else:
            g = initial_graph
            n = n - g._n
            first_node_ix = g._n
            if n <= 0:
                return g

        """
        Starting from some initial graph, we add nodes that will connect to any other node with probability 
        proportional to the indegree (or degree) of that node. Thus, diagonal-friendliness does not makes sense in
        this scenario. We also don't really care about whether the graph is directed or not since the underlying
        graph data structure will take care of this.        
        """

        # The main node generator loop
        for node_ix in range(first_node_ix, n):
            # Gather indegrees and normalize them into the target node probability.
            indegrees = g.indegrees()
            target_probabilities = indegrees / np.sum(indegrees)

            # Sample n-1 random numbers and decide whether to add a link to those nodes or not.
            r = np.random.random(target_probabilities.shape)
            edge_decision = r < target_probabilities

            # Add the node and the edges.
            g.add_node()
            for target, add_flag in enumerate(edge_decision):
                if add_flag:
                    g.add_edge(node_ix, target)

        return g

    @staticmethod
    def WattsStrogatz(n, k, beta, directed=False):
        def rewire(adj, i, j):
            possible_links = [ix for ix in np.where(adj[i] == 0) if ix != i]
            link = np.random.choice(possible_links, size=1)

            adj[i,j] = adj[j,i] = 0
            adj[i,link] = adj[link,i] = 1

        g = Network.WSRingLattice(n, k)
        adj = g._adj.copy()
        rewire_queue = []
        for i in range(n):
            for j in range(i + 1, n):
                if adj[i, j] == 1:
                    # Extract random number
                    if np.random.binomial(1, beta):
                        rewire_queue.append((i, j))

        # Apply rewiring
        for i, j in rewire_queue:
            rewire(adj, i, j)

    @staticmethod
    def ErdosRenyi(n, p, directed=False, allow_diagonals=False):
        """
        Creates a random graph.
        """

        """
        Random edge generation:
        Idea: generate an adjecency mtx with random values between 0 and 1. 
        All those over p will be added. We then filter the matrix based on whether options such as directed, 
        diagonal-friendly, etc were selected.

        To have a fair sample in the undirected case, we just mirror the upper triangle of the matrix w.r.t main diagonal
        """
        adj = np.random.random(size=(n, n))
        # Adjust for diagonal-friendliness by element-wise multiplying by (1 - identity mtx).
        # This is just swapping the 0s and 1s in the identity matrix (i.e. creating a 0-diagonal mtx).
        if not allow_diagonals:
            adj = adj * (1 - np.eye(n))

        # Adjust for directedness by mirroring the upper triangle of the matrix.
        # np.tril(n, -1) will generate indexes for the upper triangle. We then just copy the lower triangle of
        # the transpose to our original mtx.
        if not directed:
            lower_triangle_ixs = np.tril_indices(n, -1)
            adj[lower_triangle_ixs] = adj.T[lower_triangle_ixs]

        # Filter the mtx based on the probability of the edges.
        adj = (adj < p).astype(dtype=np.int8)

        # We can now safely set the new adjacency matrix.
        g = Network(directed)
        g.set_adj(adj)
        return g

    @staticmethod
    def WSRingLattice(n, k):
        """
        Constructs an undirected ring lattice with n nodes, each connected to k neighbors, k/2 on each side.
        k should be even.
        """

        adj = np.zeros((n, n))
        for i in range(n):
            for j in range(n):
                if 0 < abs(i - j) % (n - 1 - (k // 2)) <= k // 2:
                    adj[i, j] = 1

        g = Network(directed=True)
        g.set_adj(adj)
        return g

    @staticmethod
    def FullyConnected(n, directed=False, allow_diagonals=False):
        adj = np.ones(n) * (1 - np.eye(n)) + (allow_diagonals * np.eye(n))

        g = Network(directed)
        g.set_adj(adj)
        return g

    def __init__(self, directed=False):
        self._eta = 100
        self._adj = np.zeros(shape=(self._eta, self._eta), dtype=np.int8)
        self._n = 0
        self._directed = directed

    def set_adj(self, adj):
        self._adj = adj
        self._n = adj.shape[0]

    # region Analytics
    def degrees(self, kind='any'):
        if kind == 'any' or kind == 'in':
            return np.sum(self._adj, axis=0)
        elif kind == 'out':
            return np.sum(self._adj, axis=1)
        else:
            raise ValueError('Wrong degree kind requested: {}'.format(kind))

    def indegrees(self):
        return self.degrees(kind='in')

    def degree_histogram(self, kind='any'):
        degrees = self.degrees(kind)
        return np.histogram(degrees, bins=max(degrees))

    # endregion

    # region Graph manipulation
    def batch_add_node(self, n):
        for _ in range(n):
            self.add_node()

    def add_node(self):
        self._n += 1
        if self._n >= self._node_count():
            self._expand()
        return self._n - 1

    def add_edge(self, i, j):
        if (i >= self._n) or (j >= self._n) or (i < 0) or (j < 0):
            raise ValueError('Cannot add edge ({}, {})'.format(i, j))
        else:
            self._adj[i, j] = 1
            if not self._directed:
                self._adj[j, i] = 1

    # endregion

    def to_networkx(self):
        sanitized_adj = self._adj[:self._n, :self._n]

        from networkx import from_numpy_matrix
        return from_numpy_matrix(sanitized_adj)

    # region Protected
    def _node_count(self):
        return self._adj.shape[0]

    def _expand(self):
        self._adj = np.pad(self._adj, ((0, self._eta), (0, self._eta)), mode='constant', constant_values=0)
    # endregion
