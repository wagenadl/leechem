import csv

import numpy as np
import networkx as nx

from sbemdb import SBEMDB
from cleandb import clean_db


np.random.seed(0)


def distance(x1, y1, z1, x2, y2, z2):
    """Euclidean distance in 3-dimensional space between points (x1, y1, z1) and (x2, y2, z3)."""
    return np.sqrt(np.power(x1 - x2, 2) + np.power(y1 - y2, 2) + np.power(z1 - z2, 2))


def average(data):
    """Average value of iterable data structure."""
    return sum(data) / len(data)


def nearest_neighbor_distances(nodes, distance_function):
    """Calculates the distance to the nearest neighbor from the given list.

    For each node in a given list iterates through the rest of nodes,
    and finds distance to the nearest one, based on given distance function.

    Returns list of distances in same order as original list.
    I.e. element 2 represents distance from node 2 in list to its nearest neighbor.
    """

    try:
        assert len(nodes) > 1  # Nearest neighbors only makes sense for more than one node.
    except AssertionError:
        print("Actual length:", len(nodes))
        raise

    distances = list()
    for i, node_1 in enumerate(nodes):
        distances_from_node_1 = list()

        for j, node_2 in enumerate(nodes):
            if i == j:
                continue
            d = distance_function(node_1, node_2)
            distances_from_node_1.append(d)

        distances.append(min(distances_from_node_1))

    return distances


# Global variable used by functions below.
DB = clean_db(SBEMDB())


def get_node_coords():
    """Gets node coordicates of tree 444.

    DEPENDS on global 'DB'!.
    """

    x, y, z, nids = DB.nodexyz(f'tid=444')
    # nodes' coordinates map
    node_coords = {}
    for i, nid in enumerate(nids):
        node_coords[nid] = [x[i], y[i], z[i]]
    return node_coords


def get_random_sample(size):
    """Gets random node sample from tree 444.

    DEPENDS on global 'DB'!.
    """

    nodes_444 = list(DB.nodeDetails('tid=444').keys())
    sampled = np.random.choice(nodes_444, size)
    return sampled


def synapses(clause=None):
    """SYNAPSES - Find position of synapses
    (xx, yy, zz, pretid, posttid, sid_ids, prenids, postnids) = SYNAPSES(clause) returns the
    coordinates and pre- and postsynaptic tree IDs of synapses.
    CLAUSE may specify something about the presynaptic tree or the
    postsynaptic tree. E.g., 'post.tid==444'.
    If CLAUSE is not given, then all synapses connections will be retrieved.
    Results XX, YY, ZZ, expressed in microns and are corrected for
    the gapshift.
    The result is the average of the position of the presynaptic node and
    of the postsynaptic node.
    Only finds synapses that have at least one presynaptic terminal and
    at least one postsynaptic terminal.
    """

    def belowgap(_zz):
        return _zz >= 328.7

    def gapshiftx():
        return -4.285

    def gapshifty():
        return 25.857

    query = 'select a.x as ax, a.y as ay, a.z as az,' \
            + ' b.x as bx, b.y as by, b.z as bz,' \
            + ' pre.tid as atid, post.tid as btid,' \
            + ' s.sid as sid, a.nid as prenid, b.nid as postnid' \
            + ' from trees as pre' \
            + ' inner join nodes as a on pre.tid==a.tid' \
            + ' inner join syncons as sca on a.nid==sca.nid' \
            + ' inner join synapses as s on sca.sid==s.sid' \
            + ' inner join syncons as scb on scb.sid==s.sid' \
            + ' inner join nodes as b on scb.nid==b.nid' \
            + ' inner join trees as post on b.tid==post.tid'

    if clause is not None and clause != '':
        query += f' where ({clause}) and a.typ==5 and b.typ==6'
    else:
        query += f' where a.typ==5 and b.typ==6'

    c = DB.db.cursor()
    c.execute(query)

    rows = c.fetchall()
    n = len(rows)
    xx = np.zeros(n)
    yy = np.zeros(n)
    zz = np.zeros(n)
    pretid = np.zeros(n, dtype=int)
    posttid = np.zeros(n, dtype=int)
    synid = np.zeros(n, dtype=int)
    prenid = np.zeros(n, dtype=int)
    postnid = np.zeros(n, dtype=int)
    for n in range(n):
        row = rows[n]
        xx[n] = (row[0] + row[3]) / 2.0
        yy[n] = (row[1] + row[4]) / 2.0
        zz[n] = (row[2] + row[5]) / 2.0
        pretid[n] = row[6]
        posttid[n] = row[7]
        synid[n] = row[8]
        prenid[n] = row[9]
        postnid[n] = row[10]
    xx = DB.pixtoum(xx)
    yy = DB.pixtoum(yy)
    zz = DB.slicetoum(zz)
    idx = belowgap(zz)
    xx[idx] += gapshiftx()
    yy[idx] += gapshifty()
    return xx, yy, zz, pretid, posttid, synid, prenid, postnid


# Global variable used by functions below.
NODE_COORDS = get_node_coords()


def euclidean_node_distance(node_1, node_2):
    """Calculates the euclidean distance between nodes by their IDs.

    DEPENDS on global 'NODE_COORDS'!

    :param node_1: ID of first node
    :param node_2: ID of second node
    """
    return distance(*(NODE_COORDS[node_1] + NODE_COORDS[node_2]))


def get_graph():
    """Gets tree 444 as networkx.DiGraph object

    DEPENDS on globals 'DB' and 'NODE_COORDS'!.
    """

    c = DB.db.cursor()
    c.execute(f'select c.nid1, c.nid2'
              f' from nodes inner join nodecons as c on nodes.nid==c.nid1 inner join nodes as b on c.nid2==b.nid'
              f' where nodes.tid=444 and b.tid=444')

    g = nx.DiGraph()
    for row in c.fetchall():
        g_nid1, g_nid2 = row[0], row[1]
        g.add_edge(g_nid1, g_nid2, distance=distance(*(NODE_COORDS[g_nid1] + NODE_COORDS[g_nid2])))

    return g


# Global variable used by functions below.
G = get_graph()


def path_node_distance(node_1, node_2):
    """Calculates the distance between nodes by their IDs along the nearest path along the graph.

    DEPENDS on global 'G'!

    :param node_1: ID of first node
    :param node_2: ID of second node
    """
    return nx.dijkstra_path_length(G, node_1, node_2, 'distance')


def analyse_trees():

    trees_to_process = (388, 393, 32)
    sample_size = 13
    output_filepath = "tree_data.csv"

    random_sample = get_random_sample(sample_size)

    path_nn_distances = nearest_neighbor_distances(random_sample, path_node_distance)
    average_path_nn_distance = average(path_nn_distances)
    euclidean_nn_distances = nearest_neighbor_distances(random_sample, euclidean_node_distance)
    average_euclidean_nn_distance = average(euclidean_nn_distances)

    tree_data = list()
    for tree in trees_to_process:
        _, _, _, _, _, _, _, postnid = synapses(f'pre.tid={tree} and post.tid=444')

        other_tree_path_nn_distances = nearest_neighbor_distances(postnid, path_node_distance)
        other_tree_average_path_nn_distance = average(other_tree_path_nn_distances)
        other_tree_euclidean_nn_distances = nearest_neighbor_distances(postnid, euclidean_node_distance)
        other_tree_average_euclidean_nn_distance = average(other_tree_euclidean_nn_distances)

        nni_path = other_tree_average_path_nn_distance / average_path_nn_distance
        nni_euclidean = other_tree_average_euclidean_nn_distance / average_euclidean_nn_distance

        tree_data.append(
            {
                "tree_id": tree,
                "nni_path": nni_path,
                "nni_euclidean": nni_euclidean,
            }
        )

    with open(output_filepath, "w", newline="") as csvfile:
        fieldnames = ["tree_id", "nni_path", "nni_euclidean"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for tree_entry in tree_data:
            writer.writerow(tree_entry)


if __name__ == "__main__":
    analyse_trees()
