# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 18:53:39 2020
@author: Shantanu Jha, Henry Liu, Will Sun
"""
import math

from operator import mul
from fractions import Fraction
from functools import reduce
from itertools import combinations

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def nCr(n, r):
    # https://stackoverflow.com/a/3027128
    return int(reduce(mul, (Fraction(n - i, i + 1) for i in range(r)), 1))


class GraphDecoder:
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction code, and then run suitable decoders.
    """

    def __init__(self, d, T):
        self.d = d
        self.T = T
        self.S = {"X": nx.Graph(), "Z": nx.Graph()}
        self.virtual = self._specify_virtual()
        self._make_syndrome_graph()

    def _specify_virtual(self):
        virtual = {}
        virtual["X"] = []
        virtual["Z"] = []
        for j in range(0, self.d, 2):
            # X virtual nodes
            virtual["X"].append((-1, -0.5, j - 0.5))
            virtual["X"].append((-1, self.d - 0.5, j + 0.5))

            # Z virtual nodes
            virtual["Z"].append((-1, j + 0.5, -0.5))
            virtual["Z"].append((-1, j - 0.5, self.d - 0.5))
        return virtual

    def _path_degeneracy(self, a, b):
        """Calculate the number of shortest error paths that link two syndrome nodes
        through both space and time.
        Args:
            a (node): Starting or ending syndrome node (degeneracy is symmetric)
            b (node): Ending or starting syndrome node (degeneracy is symmetric)
        Raises:
            nx.exception.NodeNotFound: Nodes not both in X or Z syndrome graph
        Returns:
            int: Number of degenerate shortest paths matching this syndrome pair
        """
        # Check which subgraph node is on. If x + y is even => Z, else X.
        a_sum, b_sum = a[1] + a[2], b[1] + b[2]
        if a_sum % 2 == 0 and b_sum % 2 == 0:
            subgraph = self.S["Z"]
        elif a_sum % 2 == 1 and b_sum % 2 == 1:
            subgraph = self.S["X"]
        else:
            raise nx.exception.NodeNotFound("Nodes not both in X or Z syndrome graph.")

        return len(list(nx.all_shortest_paths(subgraph, a, b, weight="distance")))

    def _make_syndrome_graph(self):
        """
        This method injects all possible Pauli errors into the circuit for
        ``code``.
        This is done by examining the qubits used in each gate of the
        circuit for a stored logical 0. A graph is then created with a node
        for each non-trivial syndrome element, and an edge between all such
        elements that can be created by the same error.
        """

        for t in range(0, self.T):
            start_node_x = (0.5, 0.5)
            self.S["X"].add_node(
                (t,) + start_node_x,
                virtual=0,
                pos=(start_node_x[1], start_node_x[0]),
                pos_3D=(start_node_x[1], -start_node_x[0], t),
                time=t,
            )
            start_node_z = (0.5, 1.5)
            self.S["Z"].add_node(
                (t,) + start_node_z,
                virtual=0,
                pos=(start_node_z[1], start_node_z[0]),
                pos_3D=(start_node_z[1], -start_node_z[0], t),
                time=t,
            )
            self.populate_syndrome_graph((t,) + start_node_x, t, [], "X", edge_weight=1)
            self.populate_syndrome_graph((t,) + start_node_z, t, [], "Z", edge_weight=1)

        for error_key in ["X", "Z"]:
            syndrome_nodes_t0 = [
                x for x, y in self.S[error_key].nodes(data=True) if y["time"] == 0
            ]
            for node in syndrome_nodes_t0:
                space_label = (node[1], node[2])
                for t in range(0, self.T - 1):
                    self.S[error_key].add_edge(
                        (t,) + space_label, (t + 1,) + space_label, distance=1
                    )

    def populate_syndrome_graph(
        self, current_node, t, visited_nodes, error_key, edge_weight=1
    ):
        visited_nodes.append(current_node)
        neighbors = []
        i = current_node[1]
        j = current_node[2]
        neighbors.append((i - 1, j - 1))
        neighbors.append((i + 1, j - 1))
        neighbors.append((i - 1, j + 1))
        neighbors.append((i + 1, j + 1))

        normal_neighbors = [
            n
            for n in neighbors
            if self.valid_syndrome(n, error_key)
            and (t, n[0], n[1]) not in visited_nodes
        ]
        virtual_neighbors = [
            n
            for n in neighbors
            if (-1, n[0], n[1]) in self.virtual[error_key]
            and (-1, n[0], n[1]) not in visited_nodes
        ]
        #         print(error_key)
        #         print("Current: " + str(current_node))
        #         print("Normal: " + str(normal_neighbors))
        #         print("Virtual: " + str(virtual_neighbors))

        if not normal_neighbors and not virtual_neighbors:
            return

        for target in normal_neighbors:
            target_node = (t,) + target
            if not self.S[error_key].has_node(target_node):
                self.S[error_key].add_node(
                    target_node, virtual=0, pos=(target[1], target[0]), pos_3D=(target[1], -target[0], t), time=t
                )
            self.S[error_key].add_edge(current_node, target_node, distance=edge_weight)

        for target in virtual_neighbors:
            target_node = (-1,) + target
            if not self.S[error_key].has_node(target_node):
                self.S[error_key].add_node(
                    target_node, virtual=1, pos=(target[1], target[0]), pos_3D=(target[1], -target[0], -1), time=-1
                )
            self.S[error_key].add_edge(current_node, target_node, distance=edge_weight)

        for target in normal_neighbors:
            self.populate_syndrome_graph(
                (t,) + target, t, visited_nodes, error_key, edge_weight=1
            )

        for target in virtual_neighbors:
            self.populate_syndrome_graph(
                (-1,) + target, t, visited_nodes, error_key, edge_weight=1
            )

        return

    def valid_syndrome(self, node, error_key):
        i = node[0]
        j = node[1]
        if error_key == "X":
            if i > 0 and i < self.d - 1 and j < self.d and j > -1:
                return True
            else:
                return False
        elif error_key == "Z":
            if j > 0 and j < self.d - 1 and i < self.d and i > -1:
                return True
            else:
                return False

    def make_error_graph(self, nodes, error_key, err_prob=None):
        """Creates error syndrome subgraph from list of syndrome nodes. The output of
        this function is a graph that's ready for MWPM.
        If err_prob is specified, we adjust the shortest distance between syndrome
        nodes by the degeneracy of the error path.
        Args:
            nodes ([(t, x, y),]): List of changes of syndrome nodes in time.
            error_key (char): Which X/Z syndrome subgraph these nodes are from.
            err_prob (float, optional): Probability of IID data qubit X/Z flip. Defaults to None.
        Returns:
            nx.Graph: Nodes are syndromes, edges are proxy for error probabilities
        """
        virtual_dict = dict(self.S[error_key].nodes(data="virtual"))
        time_dict = dict(self.S[error_key].nodes(data="time"))
        error_graph = nx.Graph()
        nodes += self.virtual[error_key]

        for source, target in combinations(nodes, 2):
            for node in [source, target]:
                if not error_graph.has_node(node):
                    # Add all nodes to graph. We do this before adding edges so
                    # that attributes can be defined
                    error_graph.add_node(
                        node,
                        virtual=virtual_dict[node],
                        pos=(node[2], node[1]),
                        pos_3D=(node[2], -node[1], time_dict[node]),
                        time=time_dict[node],
                    )
            # Distance is proportional to the probability of this error chain, so
            # finding the maximum-weight perfect matching of the whole graph gives
            # the most likely sequence of errors that led to these syndromes.
            distance = int(
                nx.shortest_path_length(
                    self.S[error_key], source, target, weight="distance"
                )
            )

            # If err_prob is specified, we also account for path degeneracies, as:
            # ln(degeneracy) + distance * log(p / 1 - p)
            if err_prob:
                distance *= math.log(err_prob) - math.log1p(-err_prob)
                distance += math.log(self._path_degeneracy(source, target))
            else:  # Otherwise we can just assume that the log err_prob part is neg
                distance = -distance
            error_graph.add_edge(source, target, weight=distance)
        return error_graph

    def matching_graph(self, error_graph, error_key):
        time_dict = dict(self.S[error_key].nodes(data="time"))
        subgraph = nx.Graph()
        syndrome_nodes = [
            x for x, y in error_graph.nodes(data=True) if y["virtual"] == 0
        ]

        virtual_nodes = [
            x for x, y in error_graph.nodes(data=True) if y["virtual"] == 1
        ]

        for source, target in combinations(syndrome_nodes, 2):
            for node in [source, target]:
                if not subgraph.has_node(node):
                    subgraph.add_node(
                        node,
                        virtual=0,
                        pos=(node[2], node[1]),
                        pos_3D=(node[2], -node[1], time_dict[node]),
                        time=time_dict[node],
                    )
            subgraph.add_edge(
                source, target, weight=error_graph[source][target]["weight"]
            )
            
        for source in syndrome_nodes:
            potential_virtual = {}
            for target in virtual_nodes:
                potential_virtual[target] = error_graph[source][target]["weight"]   
            nearest_virtual = max(potential_virtual, key=potential_virtual.get)
            paired_virtual = nearest_virtual + source
            subgraph.add_node(
                paired_virtual, virtual=1, pos=(nearest_virtual[2], nearest_virtual[1]), pos_3D=(nearest_virtual[2], -nearest_virtual[1], nearest_virtual[0]), time=-1
            )
            subgraph.add_edge(
                source, paired_virtual, weight=potential_virtual[nearest_virtual]
            )
        
        paired_virtual_nodes = [
            x for x, y in subgraph.nodes(data=True) if y["virtual"] == 1
        ]

        for source in paired_virtual_nodes:
            for target in paired_virtual_nodes:
                subgraph.add_edge(source, target, weight=0)
        return subgraph

    def matching(self, matching_graph, error_key):
        matches = nx.max_weight_matching(matching_graph, maxcardinality=True)
        return matches
    
    def graph_3D(self, G, edge_label, angle=[-116, 22]):
        # Get node positions
        pos_3D = nx.get_node_attributes(G, 'pos_3D')
        
        # Define color range based on time
        colors = {x : plt.cm.plasma((y['time']+1)/self.T) for x,y in G.nodes(data=True)}
        # 3D network plot

        with plt.style.context(('ggplot')):
            
            fig = plt.figure(figsize=(10,7))
            ax = Axes3D(fig)
            
            # Loop on the nodes and look up in pos dictionary to extract the x,y,z coordinates of each node
            for node in G.nodes():
                xi, yi, zi = pos_3D[node]
                
                # Scatter plot
                ax.scatter(xi, yi, zi, color=colors[node], s=40*(1+G.degree(node)), edgecolors='k', alpha=0.7)
                
                # Node position label
                ax.text(xi, yi, zi, node, fontsize=14)
            
            # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
            # Those two points are the extrema of the line to be plotted
            for src, tgt in G.edges():
                x_1, y_1, z_1 = pos_3D[src]
                x_2, y_2, z_2 = pos_3D[tgt]
                
                x_line = np.array((x_1, x_2))
                y_line = np.array((y_1, y_2))
                z_line = np.array((z_1, z_2))                
            
            # Plot the connecting lines
                ax.plot(x_line, y_line, z_line, color='black', alpha=0.5)
            
            # Plot edge weight labels at midpoints
                x_mid = (x_1 + x_2)/2
                y_mid = (y_1 + y_2)/2
                z_mid = (z_1 + z_2)/2
                w = round(G[src][tgt][edge_label], 2)
                ax.text(x_mid, y_mid, z_mid, w)
        
        # Set the initial view
        ax.view_init(angle[1], angle[0])
        # Hide the axes
        ax.set_axis_off()
        
        plt.show()
             
        return
    
decoder = GraphDecoder(5,1)
G = decoder.S['Z']
decoder.graph_3D(G, 'distance')
a = [n for n in G.nodes if n == (0, 1.5, 0.5)][0]
b = [n for n in G.nodes if n == (-1,3.5,4.5)][0]
print(decoder._path_degeneracy(a, b))

decoder = GraphDecoder(5,2)
G = decoder.S['Z']
decoder.graph_3D(G,'distance')
a = [n for n in G.nodes if n == (0, 0.5, 1.5)][0]
b = [n for n in G.nodes if n == (1, 1.5, 2.5)][0]
print(decoder._path_degeneracy(a, b))

decoder = GraphDecoder(5,2)
G = decoder.S['X']
decoder.graph_3D(G,'distance')
a = [n for n in G.nodes if n == (0, 3.5, 3.5)][0]
b = [n for n in G.nodes if n == (1, 0.5, 4.5)][0]
print(decoder._path_degeneracy(a, b))

decoder = GraphDecoder(3, 1)
node_set = [(0,1.5,.5),(0,.5,1.5)]
error_graph = decoder.make_error_graph(node_set,'Z', err_prob=0.01)
decoder.graph_3D(error_graph,'weight')#note that some edges overlap

matching_graph = decoder.matching_graph(error_graph,'Z')
decoder.graph_3D(matching_graph,'weight')

g = decoder.matching(matching_graph,'Z')
for e in g:
    print(e)

decoder = GraphDecoder(9,1)
node_set = [(0,0.5,0.5),(0,1.5,1.5),(0,1.5,3.5),(0,3.5,3.5),(0,4.5,4.5),(0,4.5,6.5)]
error_graph = decoder.make_error_graph(node_set,'X', err_prob=0.01)
decoder.graph_3D(error_graph,'weight')#note that some edges overlap

matching_graph = decoder.matching_graph(error_graph,'X')
decoder.graph_3D(matching_graph,'weight')

g = decoder.matching(matching_graph,'X')
for e in g:
    print(e)

node_set = [(0,1.5,0.5),(0,0.5,1.5),(0,1.5,6.5),(0,2.5,5.5),(0,5.5,.5),(0,6.5,3.5)]
error_graph = decoder.make_error_graph(node_set,'Z', err_prob=0.01)
decoder.graph_3D(error_graph,'weight')#note that some edges overlap

matching_graph = decoder.matching_graph(error_graph,'Z')
decoder.graph_3D(matching_graph,'weight')

g = decoder.matching(matching_graph,'X')
for e in g:
    print(e)