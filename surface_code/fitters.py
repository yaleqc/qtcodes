# -*- coding: utf-8 -*-

"""
Graph decoder for surface codes
"""
import copy
import math
from itertools import combinations, product
from collections import defaultdict

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from qiskit import QuantumCircuit, execute

from .circuits import SurfaceCode

try:
    from qiskit import Aer

    HAS_AER = True
except ImportError:
    from qiskit import BasicAer

    HAS_AER = False


class GraphDecoder:
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction code, and then run suitable decoders.
    """

    def __init__(self, d, T, simulation=False):
        
        self.d = d
        self.T = T
        self.virtual = self._specify_virtual()
        self.S = {"X": nx.Graph(), "Z": nx.Graph()}
        self.simulation = simulation
        if simulation:
            self.code = SurfaceCode(d, T)
            self._make_syndrome_graph_simulate()
        else:
            self._make_syndrome_graph()

    def _specify_virtual(self):
        """Define coordinates of Z and X virtual nodes. Our convention is that Z
        virtual nodes are top/bottom and X virtual nodes are left/right.
        """
        virtual = {}
        virtual["X"] = []
        virtual["Z"] = []
        for j in range(0, self.d, 2):
            # Z virtual nodes
            virtual["Z"].append((-1, -0.5, j - 0.5))  # top
            virtual["Z"].append((-1, self.d - 0.5, j + 0.5))  # bottom

            # X virtual nodes
            virtual["X"].append((-1, j + 0.5, -0.5))  # left
            virtual["X"].append((-1, j - 0.5, self.d - 0.5))  # right
        return virtual

    def _make_syndrome_graph(self):
        start_nodes = {"Z": (0.5, 0.5), "X": (0.5, 1.5)}
        for error_key in ["X", "Z"]:
            # subgraphs for each time step
            for t in range(0, self.T):
                start_node = start_nodes[error_key]
                self.S[error_key].add_node(
                    (t,) + start_node,
                    virtual=0,
                    pos=(start_node[1], -start_node[0]),
                    time=t,
                    pos_3D=(
                        start_node[1],
                        -start_node[0],
                        t,
                    ),  # y-coord is flipped for plot purposes
                )
                self.populate_syndrome_graph(
                    (t,) + start_node, t, [], error_key, edge_weight=1
                )

            # connect physical qubits in same location across subgraphs of adjacent times
            syndrome_nodes_t0 = [
                x for x, y in self.S[error_key].nodes(data=True) if y["time"] == 0
            ]
            for node in syndrome_nodes_t0:
                space_label = (node[1], node[2])
                for t in range(0, self.T - 1):
                    self.S[error_key].add_edge(
                        (t,) + space_label, (t + 1,) + space_label, distance=1
                    )

    def _make_syndrome_graph_simulate(self):
        """
        This method injects all possible Pauli errors into the circuit for
        ``code``.

        This is done by examining the qubits used in each gate of the
        circuit for a stored logical 0. A graph is then created with a node
        for each non-trivial syndrome element, and an edge between all such
        elements that can be created by the same error.
        """

        qc = self.code.circuit["0"]

        blank_qc = QuantumCircuit()
        for qreg in qc.qregs:
            blank_qc.add_register(qreg)
        for creg in qc.cregs:
            blank_qc.add_register(creg)

        error_circuit = {}
        circuit_name = {}
        depth = len(qc)
        for j in range(depth):
            qubits = qc.data[j][1]
            for qubit in qubits:
                for error in ["x", "z"]:
                    temp_qc = copy.deepcopy(blank_qc)
                    temp_qc.name = str((j, qubit, error))
                    temp_qc.data = qc.data[0:j]
                    getattr(temp_qc, error)(qubit)
                    temp_qc.data += qc.data[j : depth + 1]
                    circuit_name[(j, qubit, error)] = temp_qc.name
                    error_circuit[temp_qc.name] = temp_qc

        if HAS_AER:
            simulator = Aer.get_backend("qasm_simulator")
        else:
            simulator = BasicAer.get_backend("qasm_simulator")

        job = execute(list(error_circuit.values()), simulator)

        for j in range(depth):
            qubits = qc.data[j][1]
            for qubit in qubits:
                for error in ["x", "z"]:
                    raw_results = {}
                    raw_results["0"] = job.result().get_counts(str((j, qubit, error)))

                    results = self.code.process_results(raw_results["0"])
                    extracted_nodes = self.code.extract_nodes(results)

                    for err_key, nodes in dict(
                        zip(("X", "Z"), extracted_nodes)
                    ).items():
                        # Add virtual syndrome nodes
                        for node in self.virtual[err_key]:
                            # Visualization coords (y-coord flipped for plot)
                            pos_2D = (node[2], -node[1])
                            pos_3D = (*pos_2D, (self.T - 1) / 2)  # Plot midway in stack
                            self.S[err_key].add_node(
                                node, virtual=1, pos=pos_2D, time=-1, pos_3D=pos_3D,
                            )

                        # Add syndrome nodes
                        for node in nodes:
                            # Visualization coords (y-coord flipped for plot)
                            pos_2D = (node[2], -node[1])
                            pos_3D = (*pos_2D, node[0])
                            self.S[err_key].add_node(
                                node,
                                virtual=0,
                                pos=pos_2D,
                                time=node[0],
                                pos_3D=pos_3D,
                            )

                            # Check if any neighbors are virtual nodes
                            candidates = [
                                (-1, node[1] + di, node[2] + dj)
                                for di, dj in product((1, -1), repeat=2)
                            ]
                            for virtual in candidates:
                                if virtual in self.virtual[err_key]:
                                    self.S[err_key].add_edge(node, virtual, distance=1)

                        # Add connections
                        for source, target in combinations(nodes, 2):
                            self.S[err_key].add_edge(source, target, distance=1)

    def populate_syndrome_graph(
        self, current_node, t, visited_nodes, error_key, edge_weight=1
    ):
        """Recursive function to populate syndrome subgraph at time t with error_key X/Z. The current_node
        is connected to neighboring nodes without revisiting a node.

        Args:
            current_node ((t, x, y)): Current syndrome node to be connected with neighboring nodes.
            visited_nodes ([(t, x, y),]): List of syndrome nodes which have already been traver.
            error_key (char): Which X/Z syndrome subgraph these nodes are from.
            edge_weight (float, optional): Weight of edge between two adjacent syndrome nodes. Defaults to 1.

        Returns:
            None: function is to traverse the syndrome nodes and connect neighbors
        """
        visited_nodes.append(current_node)
        neighbors = []
        i = current_node[1]  # syndrome node x coordinate
        j = current_node[2]  # syndrome node y coordinate
        neighbors.append((i - 1, j - 1))  # up left
        neighbors.append((i + 1, j - 1))  # down left
        neighbors.append((i - 1, j + 1))  # up right
        neighbors.append((i + 1, j + 1))  # down right

        normal_neighbors = [
            n
            for n in neighbors
            if self.valid_syndrome(n, error_key)
            and (t, n[0], n[1]) not in visited_nodes
        ]  # syndrome node neighbors of current_node not already visited
        virtual_neighbors = [
            n
            for n in neighbors
            if (-1, n[0], n[1]) in self.virtual[error_key]
            and (-1, n[0], n[1]) not in visited_nodes
        ]  # virtual node neighbors of current_node not already visited

        # no neighbors to add edges
        if not normal_neighbors and not virtual_neighbors:
            return

        # add normal/non-virtual neighbors
        for target in normal_neighbors:
            target_node = (
                t,
            ) + target  # target_node has time t with x and y coordinates from target
            if not self.S[error_key].has_node(target_node):
                self.S[error_key].add_node(
                    target_node,
                    virtual=0,
                    pos=(target[1], -target[0]),
                    time=t,
                    pos_3D=(target[1], -target[0], t),
                )  # add target_node to syndrome subgraph if it doesn't already exist
            self.S[error_key].add_edge(
                current_node, target_node, distance=edge_weight
            )  # add edge between current_node and target_node

        # add virtual neighbors
        for target in virtual_neighbors:
            target_node = (
                -1,
            ) + target  # virtual target_node has time -1 with x and y coordinates from target
            if not self.S[error_key].has_node(target_node):
                self.S[error_key].add_node(
                    target_node,
                    virtual=1,
                    pos=(target[1], -target[0]),
                    time=-1,
                    pos_3D=(target[1], -target[0], (self.T - 1) / 2),
                )  # add virtual target_node to syndrome subgraph with z coordinate (T-1)/2 for nice plotting, if it doesn't already exist
            self.S[error_key].add_edge(
                current_node, target_node, distance=edge_weight
            )  # add edge between current_node and virtual target_node

        # recursively traverse normal neighbors
        for target in normal_neighbors:
            self.populate_syndrome_graph(
                (t,) + target, t, visited_nodes, error_key, edge_weight=1
            )

        # recursively traverse virtual neighbors
        for target in virtual_neighbors:
            self.populate_syndrome_graph(
                (-1,) + target, t, visited_nodes, error_key, edge_weight=1
            )

    def valid_syndrome(self, node, error_key):
        """Checks whether a node is a syndrome node under our error_key, which is either X or Z.

        Args:
            node ((t, x, y)): Node in graph.
            error_key (char): Which X/Z syndrome subgraph these nodes are from.

        Returns:
            Boolean T/F: whether node is a syndrome node
        """
        i = node[0]
        j = node[1]
        if error_key == "Z":
            if i > 0 and i < self.d - 1 and j < self.d and j > -1:
                return True
            else:
                return False
        elif error_key == "X":
            if j > 0 and j < self.d - 1 and i < self.d and i > -1:
                return True
            else:
                return False

    def make_error_graph(self, nodes, error_key, err_prob=None):
        """Creates error syndrome subgraph from list of syndrome nodes. The output of
        this function is a graph that's ready for minimum weight perfect matching (MWPM).

        If err_prob is specified, we adjust the shortest distance between syndrome
        nodes by the degeneracy of the error path.

        Args:
            nodes ([(t, x, y),]): List of changes of syndrome nodes in time.
            error_key (char): Which X/Z syndrome subgraph these nodes are from.
            err_prob (float, optional): Probability of IID data qubit X/Z flip. Defaults to None.

        Returns:
            nx.Graph: Nodes are syndromes, edges are proxy for error probabilities
        """
        paths = {}
        virtual_dict = nx.get_node_attributes(self.S[error_key], "virtual")
        time_dict = nx.get_node_attributes(self.S[error_key], "time")
        error_graph = nx.Graph()
        nodes += self.virtual[error_key]

        for node in nodes:
            if not error_graph.has_node(node):
                error_graph.add_node(
                    node,
                    virtual=virtual_dict[node],
                    pos=(node[2], -node[1]),
                    time=time_dict[node],
                    pos_3D=(node[2], -node[1], time_dict[node]),
                )

        for source, target in combinations(nodes, 2):
            # Distance is proportional to the probability of this error chain, so
            # finding the maximum-weight perfect matching of the whole graph gives
            # the most likely sequence of errors that led to these syndromes.
            distance = int(
                nx.shortest_path_length(
                    self.S[error_key], source, target, weight="distance"
                )
            )

            # If err_prob is specified, we also account for path degeneracies
            deg, path = self._path_degeneracy(source, target, error_key)
            paths[(source, target)] = path
            if err_prob:
                distance = distance - math.log(deg)/(math.log1p(-err_prob) - math.log(err_prob))
            distance = -distance    
            error_graph.add_edge(source, target, weight=distance)


        if self.simulation: #paths incorrect for simulated syndrome graph
            return error_graph

        return error_graph, paths

    def analytic_paths(self, matches, error_key):
        analytic_decoder = GraphDecoder(self.d,self.T)
        paths = {}
        for (source,target) in matches:
            _, path = analytic_decoder._path_degeneracy(source[:3],target[:3], error_key)
            paths[(source[:3], target[:3])] = path
        return paths

    def _path_degeneracy(self, a, b, error_key):
        """Calculate the number of shortest error paths that link two syndrome nodes
        through both space and time.

        Args:
            a (tuple): Starting or ending syndrome node (degeneracy is symmetric)
            b (tuple): Ending or starting syndrome node (degeneracy is symmetric)

        Raises:
            nx.exception.NodeNotFound: error_key must be X or Z

        Returns:
            int: Number of degenerate shortest paths matching this syndrome pair
            [nodes,]: List of nodes for one of the shortest paths
        """
        # Check which subgraph node is on. If x + y is even => X, else Z.
        # a_sum, b_sum = a[1] + a[2], b[1] + b[2]
        if error_key == "X":
            subgraph = self.S["X"]
        elif error_key == "Z":
            subgraph = self.S["Z"]
        else:
            raise nx.exception.NodeNotFound("error_key must be X or Z")

        shortest_paths = list(nx.all_shortest_paths(subgraph, a, b, weight="distance"))
        one_path = shortest_paths[
            0
        ]  # We can pick any path to return as the error chain
        degeneracy = len(shortest_paths)

        # If either node is a virtual node, we also find degeneracies from the other
        # node to *any* nearest virtual node
        source = None
        if a[0] == -1:
            target = a
            source = b
        elif b[0] == -1:
            target = b
            source = a

        # Compute additional degeneracies to edge boundaries
        if source:
            virtual_nodes = self.virtual[error_key]
            shortest_distance = nx.shortest_path_length(
                subgraph, a, b, weight="distance"
            )
            for node in virtual_nodes:
                distance = nx.shortest_path_length(
                    subgraph, source, node, weight="distance"
                )
                if distance == shortest_distance and node != target:
                    degeneracy += len(
                        list(
                            nx.all_shortest_paths(
                                subgraph, source, node, weight="distance"
                            )
                        )
                    )
        return degeneracy, one_path

    def matching_graph(self, error_graph, error_key):
        """Return subgraph of error graph to be matched.

        Args:
            error_graph (nx.Graph): Complete error graph to be matched.
            error_key (char): Which X/Z syndrome subgraph these nodes are from.

        Returns:
            nx.Graph: Subgraph of error graph to be matched
        """
        time_dict = nx.get_node_attributes(self.S[error_key], "time")
        subgraph = nx.Graph()
        syndrome_nodes = [
            x for x, y in error_graph.nodes(data=True) if y["virtual"] == 0
        ]
        virtual_nodes = [
            x for x, y in error_graph.nodes(data=True) if y["virtual"] == 1
        ]

        # add and connect each syndrome node to subgraph
        for node in syndrome_nodes:
            if not subgraph.has_node(node):
                subgraph.add_node(
                    node,
                    virtual=0,
                    pos=(node[2], -node[1]),
                    time=time_dict[node],
                    pos_3D=(node[2], -node[1], time_dict[node]),
                )
        for source, target in combinations(syndrome_nodes, 2):
            subgraph.add_edge(
                source, target, weight=error_graph[source][target]["weight"]
            )

        # connect each syndrome node to its closest virtual node in subgraph
        for source in syndrome_nodes:
            potential_virtual = {}
            for target in virtual_nodes:
                potential_virtual[target] = error_graph[source][target]["weight"]
            nearest_virtual = max(potential_virtual, key=potential_virtual.get)
            paired_virtual = (
                nearest_virtual + source
            )  # paired_virtual (virtual, syndrome) allows for the virtual node to be matched more than once
            subgraph.add_node(
                paired_virtual,
                virtual=1,
                pos=(nearest_virtual[2], -nearest_virtual[1]),
                time=-1,
                pos_3D=(nearest_virtual[2], -nearest_virtual[1], -1),
            )  # add paired_virtual to subgraph
            subgraph.add_edge(
                source, paired_virtual, weight=potential_virtual[nearest_virtual]
            )  # add (syndrome, paired_virtual) edge to subgraph

        paired_virtual_nodes = [
            x for x, y in subgraph.nodes(data=True) if y["virtual"] == 1
        ]

        # add 0 weight between paired virtual nodes
        for source, target in combinations(paired_virtual_nodes, 2):
            subgraph.add_edge(source, target, weight=0)

        return subgraph

    def matching(self, matching_graph, error_key):
        """Return matches of minimum weight perfect matching (MWPM) on matching_graph.

        Args:
            matching_graph (nx.Graph): Graph to run MWPM.
            error_key (char): Which X/Z syndrome subgraph these nodes are from.

        Returns:
            [(node, node),]: List of matchings found from MWPM
        """
        matches = nx.max_weight_matching(matching_graph, maxcardinality=True)
        filtered_matches = [
            (source, target)
            for (source, target) in matches
            if not (len(source) > 3 and len(target) > 3)
        ]  # remove 0 weighted matched edges between virtual syndrome nodes
        return filtered_matches

    def calculate_qubit_flips(self, matches, paths, error_key):
        physical_qubit_flips = {}
        for (source, target) in matches:
            # Trim "paired virtual" nodes to nearest virtual node
            if len(source) > 3:
                source = source[:3]
            if len(target) > 3:
                target = target[:3]

            # Paths dict is encoded in one direction, check other if not found
            if (source, target) not in paths:
                source, target = (target, source)

            path = paths[(source, target)]  # This is an arbitrary shortest error path
            for i in range(0, len(path) - 1):
                start = path[i]
                end = path[i + 1]
                # Check if syndromes are in different physical locations
                # If they're in the same location, this is a measurement error
                if start[1:] != end[1:]:
                    time = start[0]
                    if time == -1:  # Grab time from non-virtual syndrome
                        time = end[0]
                    physical_qubit = (
                        time,
                        (start[1] + end[1]) / 2,
                        (start[2] + end[2]) / 2,
                    )

                    # Paired flips at the same time can be ignored
                    if physical_qubit in physical_qubit_flips:
                        physical_qubit_flips[physical_qubit] = (
                            physical_qubit_flips[physical_qubit] + 1
                        ) % 2
                    else:
                        physical_qubit_flips[physical_qubit] = 1

        physical_qubit_flips = {
            x: error_key for x, y in physical_qubit_flips.items() if y == 1
        }
        return physical_qubit_flips

    def net_qubit_flips(self, flips_x, flips_z):
        flipsx = {flip: "X" for flip, _ in flips_x.items() if flip not in flips_z}
        flipsz = {flip: "Z" for flip, _ in flips_z.items() if flip not in flips_x}
        flipsy = {flip: "Y" for flip, _ in flips_x.items() if flip in flips_z}
        flips = {**flipsx, **flipsy, **flipsz}

        individual_flips = defaultdict(dict)

        for flip, error_key in flips.items():
            individual_flips[flip[1:]][flip[0]] = error_key

        paulis = {
            "X": np.array([[0, 1], [1, 0]]),
            "Y": np.array([[0, -1j], [1j, 0]]),
            "Z": np.array([[1, 0], [0, -1]]),
            "I": np.array([[1, 0], [0, 1]]),
        }

        physical_qubit_flips = {}
        for qubit_loc, flip_record in individual_flips.items():
            net_error = paulis["I"]
            # print("Physical Qubit: " + str(qubit_loc))
            for time, error in sorted(flip_record.items(), key=lambda item: item[0]):
                # print("Error: " + error + " at time: " + str(time))
                net_error = net_error.dot(paulis[error])
            physical_qubit_flips[qubit_loc] = net_error

        physical_qubit_flips = {x:y for x,y in physical_qubit_flips.items() if not np.array_equal(y,paulis["I"])}
        return physical_qubit_flips

    def graph_2D(self, G, edge_label):
        pos = nx.get_node_attributes(G, "pos")
        nx.draw_networkx(G, pos)
        labels = nx.get_edge_attributes(G, edge_label)
        labels = {x: round(y, 3) for (x, y) in labels.items()}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
        plt.show()

    def graph_3D(self, G, edge_label, angle=[-116, 22]):
        """Plots a graph with edge labels in 3D.

        Args:
            G (nx.Graph): Graph to plot in 3D.
            edge_label (float): Edge label to display; either distance or weight.
            angle ([float, float]): Initial 3D angle view. Defaults to [-116, 22]

        Returns:
            None: Plot is displayed in plt.show()
        """
        # Get node 3D positions
        pos_3D = nx.get_node_attributes(G, "pos_3D")

        # Define color range based on time
        colors = {
            x: plt.cm.plasma((y["time"] + 1) / self.T) for x, y in G.nodes(data=True)
        }

        # 3D network plot
        with plt.style.context(("ggplot")):

            fig = plt.figure(figsize=(20, 14))
            ax = Axes3D(fig)

            # Loop on the nodes and look up in pos dictionary to extract the x,y,z coordinates of each node
            for node in G.nodes():
                xi, yi, zi = pos_3D[node]

                # Scatter plot
                ax.scatter(
                    xi,
                    yi,
                    zi,
                    color=colors[node],
                    s=120 * (1 + G.degree(node)),
                    edgecolors="k",
                    alpha=0.7,
                )

                # Label node position
                ax.text(xi, yi, zi, node, fontsize=20)

            # Loop on the edges to get the x,y,z, coordinates of the connected nodes
            # Those two points are the extrema of the line to be plotted
            for src, tgt in G.edges():
                x_1, y_1, z_1 = pos_3D[src]
                x_2, y_2, z_2 = pos_3D[tgt]

                x_line = np.array((x_1, x_2))
                y_line = np.array((y_1, y_2))
                z_line = np.array((z_1, z_2))

                # Plot the connecting lines
                ax.plot(x_line, y_line, z_line, color="black", alpha=0.5)

                # Label edges at midpoints
                x_mid = (x_1 + x_2) / 2
                y_mid = (y_1 + y_2) / 2
                z_mid = (z_1 + z_2) / 2
                label = round(G[src][tgt][edge_label], 2)
                ax.text(x_mid, y_mid, z_mid, label, fontsize=14)

        # Set the initial view
        ax.view_init(angle[1], angle[0])

        # Hide the axes
        ax.set_axis_off()

        # Get rid of colored axes planes
        # First remove fill
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False

        # Now set color to white (or whatever is "invisible")
        ax.xaxis.pane.set_edgecolor("w")
        ax.yaxis.pane.set_edgecolor("w")
        ax.zaxis.pane.set_edgecolor("w")

        plt.show()
