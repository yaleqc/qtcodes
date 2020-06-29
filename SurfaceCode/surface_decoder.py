# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 18:53:39 2020

@author: Shantanu Jha
"""
from operator import mul
from fractions import Fraction
from functools import reduce

import networkx as nx


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

        return len(list(nx.all_shortest_paths(subgraph, a, b)))

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
                time=t,
            )
            start_node_z = (0.5, 1.5)
            self.S["Z"].add_node(
                (t,) + start_node_z,
                virtual=0,
                pos=(start_node_z[1], start_node_z[0]),
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
                    target_node, virtual=0, pos=(target[1], target[0]), time=t
                )
            self.S[error_key].add_edge(current_node, target_node, distance=edge_weight)

        for target in virtual_neighbors:
            target_node = (-1,) + target
            if not self.S[error_key].has_node(target_node):
                self.S[error_key].add_node(
                    target_node, virtual=1, pos=(target[1], target[0]), time=-1
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

    def make_error_graph(self, nodes, error_key):
        virtual_dict = dict(self.S[error_key].nodes(data="virtual"))
        time_dict = dict(self.S[error_key].nodes(data="time"))
        error_graph = nx.Graph()
        nodes += self.virtual[error_key]
        for source in nodes:
            for target in nodes:
                if source != target:
                    for node in [source, target]:
                        if not error_graph.has_node(node):
                            error_graph.add_node(
                                node,
                                virtual=virtual_dict[node],
                                pos=(node[2], node[1]),
                                time=time_dict[node],
                            )
                    distance = int(
                        nx.shortest_path_length(
                            self.S[error_key], source, target, weight="distance"
                        )
                    )
                    error_graph.add_edge(source, target, weight=-distance)
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

        for source in syndrome_nodes:
            for target in syndrome_nodes:
                if source != target:
                    for node in [source, target]:
                        if not subgraph.has_node(node):
                            subgraph.add_node(
                                node,
                                virtual=0,
                                pos=(node[2], node[1]),
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
                paired_virtual, virtual=1, pos=(nearest_virtual[2], nearest_virtual[1])
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
