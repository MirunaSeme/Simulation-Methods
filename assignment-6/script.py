import os

from bokeh.models.graphs import from_networkx

from network import Network
from plot import plot_networkx

if __name__ == '__main__':
    g = Network.AlbertBarabasi(100)
    plot_networkx(g.to_networkx())