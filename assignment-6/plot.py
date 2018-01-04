import networkx as nx

from bokeh.io import show, output_file
from bokeh.plotting import figure
from bokeh.models.graphs import from_networkx

import matplotlib.pyplot as plt

# TODO: switch to networkx plotting. Try and do size/color encoding of degree/indegree/outdegree, etc.
# TODO: try to make graph interactive



def plot_networkx_through_bokeh(g, graph_name='MyGraph'):
    plot = figure(title=graph_name,
                  x_range=(-1.1, 1.1),
                  y_range=(-1.1, 1.1),
                  plot_width=1000,
                  plot_height=1000)

    graph = from_networkx(g,
                          nx.spring_layout,
                          scale=3,
                          center=(0, 0),
                          iterations=200)

    plot.renderers.append(graph)

    output_file("networkx_graph.html")
    show(plot)


def plot_networkx(g):
    nodes, degrees = zip(*nx.degree(g))
    nx.draw(g, nodelist=nodes, node_size=[v * 100 for v in degrees])
    plt.show()