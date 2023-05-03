import math

import networkx


class EuclideanNodeDistance:

    @staticmethod
    def distance(node_1, node_2):

        return math.sqrt(
            (node_1.x - node_2.x) ** 2
            + (node_1.y - node_2.y) ** 2
            + (node_1.z - node_2.z) ** 2
        )


class GraphNodeDistance:
    def __init__(self, graph):
        self._graph = graph

    def distance(self, node_1, node_2):
        return networkx.dijkstra_path_length(self._graph, node_1, node_2, 'distance')
