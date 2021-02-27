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
from typing import Tuple, List, Dict, Optional, TypeVar, Union, cast
from .base import TopologicalGraphDecoder
from .visualization import VisualizationMixin
from .mwpm import MWPMDecodingMixin

from qiskit import QuantumCircuit, execute

try:
    from qiskit import Aer

    HAS_AER = True
except ImportError:
    from qiskit import BasicAer

    HAS_AER = False


TQubit = Tuple[float, float, float]  # (time,row,column) ==> (t,i,j)
TQubitLoc = Tuple[float, float]  # (row,column) ==> (i,j)


class XXZZGraphDecoderBase(TopologicalGraphDecoder[TQubit]):
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction code, and then run suitable decoders.
    """

    def __init__(self, code_params: Dict):
        super().__init__(code_params)
        if "d" not in self.code_params or "T" not in self.code_params:
            raise ValueError("Please include d and T in code_params.")

        self.virtual = self._specify_virtual()
        self.S["X"] = nx.Graph()
        self.S["Z"] = nx.Graph()
        self._make_syndrome_graph()

        self.paulis = {
            "X": np.array([[0, 1], [1, 0]]),
            "Y": np.array([[0, -1j], [1j, 0]]),
            "Z": np.array([[1, 0], [0, -1]]),
            "I": np.array([[1, 0], [0, 1]]),
        }

    def _specify_virtual(self) -> Dict[str, List[TQubit]]:
        """
        Define coordinates of Z and X virtual syndrome nodes. Our convention is that Z
        virtual syndrome nodes are top/bottom and X virtual nodes are left/right.
        Args:
        Returns:
            virtual (dictionary): where virtual["X"] holds a list of tuples specifying virtual X syndrome nodes
            and equivalently for virtual["Z"]
        """
        virtual: Dict[str, List[TQubit]] = {}
        virtual["X"] = []
        virtual["Z"] = []
        for j in range(0, self.code_params["d"], 2):
            # Z virtual nodes
            virtual["Z"].append((-1, -0.5, j - 0.5))  # top
            virtual["Z"].append((-1, self.code_params["d"] - 0.5, j + 0.5))  # bottom

            # X virtual nodes
            virtual["X"].append((-1, j + 0.5, -0.5))  # left
            virtual["X"].append((-1, j - 0.5, self.code_params["d"] - 0.5))  # right
        return virtual

    def _make_syndrome_graph(self) -> None:
        """
        Populates self.S["X"] and self.S["Z"] syndrome nx.Graph()'s with nodes specified by time and position. 
        Args:
        Returns:
        """
        start_nodes = {"Z": (0.5, 0.5), "X": (0.5, 1.5)}
        for syndrome_graph_key in ["X", "Z"]:
            # subgraphs for each time step
            for t in range(0, self.code_params["T"]):
                start_node = start_nodes[syndrome_graph_key]
                self.S[syndrome_graph_key].add_node(
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
                self._populate_syndrome_graph(
                    (t,) + start_node, t, [], syndrome_graph_key, edge_weight=1
                )

            # connect physical qubits in same location across subgraphs of adjacent times
            syndrome_nodes_t0 = [
                x
                for x, y in self.S[syndrome_graph_key].nodes(data=True)
                if y["time"] == 0
            ]
            for node in syndrome_nodes_t0:
                space_label = (node[1], node[2])
                for t in range(0, self.code_params["T"] - 1):
                    self.S[syndrome_graph_key].add_edge(
                        (t,) + space_label, (t + 1,) + space_label, distance=1
                    )

    def _populate_syndrome_graph(
        self,
        current_node: TQubit,
        t: int,
        visited_nodes: List[TQubit],
        syndrome_graph_key: str,
        edge_weight: int = 1,
    ) -> None:
        """Recursive function to populate syndrome subgraph at time t with syndrome_graph_key X/Z. The current_node
        is connected to neighboring nodes without revisiting a node.

        Args:
            current_node ((t, x, y)): Current syndrome node to be connected with neighboring nodes.
            visited_nodes ([(t, x, y),]): List of syndrome nodes which have already been traver.
            syndrome_graph_key (char): Which X/Z syndrome subgraph these nodes are from.
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
            if self._valid_syndrome(n, syndrome_graph_key)
            and (t, n[0], n[1]) not in visited_nodes
        ]  # syndrome node neighbors of current_node not already visited
        virtual_neighbors = [
            n
            for n in neighbors
            if (-1, n[0], n[1]) in self.virtual[syndrome_graph_key]
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
            if not self.S[syndrome_graph_key].has_node(target_node):
                self.S[syndrome_graph_key].add_node(
                    target_node,
                    virtual=0,
                    pos=(target[1], -target[0]),
                    time=t,
                    pos_3D=(target[1], -target[0], t),
                )  # add target_node to syndrome subgraph if it doesn't already exist
            self.S[syndrome_graph_key].add_edge(
                current_node, target_node, distance=edge_weight
            )  # add edge between current_node and target_node

        # add virtual neighbors
        for target in virtual_neighbors:
            target_node = (
                -1,
            ) + target  # virtual target_node has time -1 with x and y coordinates from target
            if not self.S[syndrome_graph_key].has_node(target_node):
                self.S[syndrome_graph_key].add_node(
                    target_node,
                    virtual=1,
                    pos=(target[1], -target[0]),
                    time=-1,
                    pos_3D=(target[1], -target[0], (self.code_params["T"] - 1) / 2),
                )  # add virtual target_node to syndrome subgraph with z coordinate (T-1)/2 for nice plotting, if it doesn't already exist
            self.S[syndrome_graph_key].add_edge(
                current_node, target_node, distance=edge_weight
            )  # add edge between current_node and virtual target_node

        # recursively traverse normal neighbors
        for target in normal_neighbors:
            self._populate_syndrome_graph(
                (t,) + target, t, visited_nodes, syndrome_graph_key, edge_weight=1
            )

        # recursively traverse virtual neighbors
        for target in virtual_neighbors:
            self._populate_syndrome_graph(
                (-1,) + target, t, visited_nodes, syndrome_graph_key, edge_weight=1
            )

    def _valid_syndrome(self, node: TQubitLoc, syndrome_graph_key: str) -> bool:
        """Checks whether a node is a syndrome node under our syndrome_graph_key, which is either X or Z.

        Args:
            node ((t, x, y)): Node in graph.
            syndrome_graph_key (char): Which X/Z syndrome subgraph these nodes are from.

        Returns:
            Boolean T/F: whether node is a syndrome node
        """
        i = node[0]
        j = node[1]
        if syndrome_graph_key == "Z":
            if (
                i > 0
                and i < self.code_params["d"] - 1
                and j < self.code_params["d"]
                and j > -1
            ):
                return True
            else:
                return False
        elif syndrome_graph_key == "X":
            if (
                j > 0
                and j < self.code_params["d"] - 1
                and i < self.code_params["d"]
                and i > -1
            ):
                return True
            else:
                return False
        else:
            raise ValueError("Please enter a valid syndrome_graph_key: X or Z")

    def _make_error_graph(
        self,
        nodes: List[TQubit],
        syndrome_graph_key: str,
        err_prob: Optional[int] = None,
    ):
        """Creates error syndrome subgraph from list of syndrome nodes. The output of
        this function is a graph that's ready for minimum weight perfect matching (MWPM).

        If err_prob is specified, we adjust the shortest distance between syndrome
        nodes by the degeneracy of the error path.

        Args:
            nodes ([(t, x, y),]): List of changes of syndrome nodes in time.
            syndrome_graph_key (char): Which X/Z syndrome subgraph these nodes are from.
            err_prob (float, optional): Probability of IID data qubit X/Z flip. Defaults to None.

        Returns:
            nx.Graph: Nodes are syndromes, edges are proxy for error probabilities
        """
        paths = {}
        virtual_dict = nx.get_node_attributes(self.S[syndrome_graph_key], "virtual")
        time_dict = nx.get_node_attributes(self.S[syndrome_graph_key], "time")
        error_graph = nx.Graph()
        nodes += self.virtual[syndrome_graph_key]

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
            distance = float(
                nx.shortest_path_length(
                    self.S[syndrome_graph_key], source, target, weight="distance"
                )
            )

            # If err_prob is specified, we also account for path degeneracies
            deg, path = self._path_degeneracy(source, target, syndrome_graph_key)
            paths[(source, target)] = path
            if err_prob:
                distance = distance - math.log(deg) / (
                    math.log1p(-err_prob) - math.log(err_prob)
                )
            distance = -distance
            error_graph.add_edge(source, target, weight=distance)

        return error_graph, paths

    def _path_degeneracy(self, a: TQubit, b: TQubit, syndrome_graph_key: str):
        """Calculate the number of shortest error paths that link two syndrome nodes
        through both space and time.

        Args:
            a (tuple): Starting or ending syndrome node (degeneracy is symmetric)
            b (tuple): Ending or starting syndrome node (degeneracy is symmetric)

        Raises:
            nx.exception.NodeNotFound: syndrome_graph_key must be X or Z

        Returns:
            int: Number of degenerate shortest paths matching this syndrome pair
            [nodes,]: List of nodes for one of the shortest paths
        """
        # Check which subgraph node is on. If x + y is even => X, else Z.
        # a_sum, b_sum = a[1] + a[2], b[1] + b[2]
        if syndrome_graph_key == "X":
            subgraph = self.S["X"]
        elif syndrome_graph_key == "Z":
            subgraph = self.S["Z"]
        else:
            raise nx.exception.NodeNotFound("syndrome_graph_key must be X or Z")

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
            virtual_nodes = self.virtual[syndrome_graph_key]
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

    def _make_matching_graph(
        self, error_graph: nx.Graph, syndrome_graph_key: str
    ) -> nx.Graph:
        """Return subgraph of error graph to be matched.

        Args:
            error_graph (nx.Graph): Complete error graph to be matched.
            syndrome_graph_key (char): Which X/Z syndrome subgraph these nodes are from.

        Returns:
            nx.Graph: Subgraph of error graph to be matched
        """
        time_dict = nx.get_node_attributes(self.S[syndrome_graph_key], "time")
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
            potential_virtual: Dict[TQubit, float] = {}
            for target in virtual_nodes:
                potential_virtual[target] = error_graph[source][target]["weight"]
            nearest_virtual = max(potential_virtual, key=lambda k: potential_virtual[k])
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

    def _run_mwpm(
        self, matching_graph: nx.Graph, syndrome_graph_key: str,
    ) -> List[Tuple[TQubit, TQubit]]:
        """Return matches of minimum weight perfect matching (MWPM) on matching_graph.

        Args:
            matching_graph (nx.Graph): Graph to run MWPM.
            syndrome_graph_key (char): Which X/Z syndrome subgraph these nodes are from.

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

    def _get_matched_graph(
        self, matching_graph: nx.Graph, filtered_matches: List[Tuple[TQubit, TQubit]],
    ) -> nx.Graph:
        matched_graph = matching_graph.copy()
        for u, v, _ in matching_graph.edges(data=True):
            if (u, v) not in filtered_matches and (v, u) not in filtered_matches:
                matched_graph.remove_edge(u, v)
        matched_graph.remove_nodes_from(list(nx.isolates(matched_graph)))
        return matched_graph

    def _calculate_qubit_flips(
        self,
        matches: List[Tuple[TQubit, TQubit]],
        paths: Dict[Tuple[TQubit, TQubit], List[TQubit]],
        syndrome_graph_key,
    ) -> List[TQubit]:
        """Return all data qubits with odd syndrome_graph_key flips

        Args:
            matches ([(node,node),]): MWPM chosen node pairs
            paths ({(node,node):[node,]}): shortest path of nodes for each node pair
            syndrome_graph_key (char): Which X/Z syndrome subgraph these nodes are from.

        Returns:
            [node,]: List of data qubits with odd syndrome_graph_key flips
        """
        data_qubit_flips: Dict[TQubit, int] = {}
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
                    data_qubit = (
                        time,
                        (start[1] + end[1]) / 2,
                        (start[2] + end[2]) / 2,
                    )

                    # Paired flips at the same time can be ignored
                    if data_qubit in data_qubit_flips:
                        data_qubit_flips[data_qubit] = (
                            data_qubit_flips[data_qubit] + 1
                        ) % 2
                    else:
                        data_qubit_flips[data_qubit] = 1

        data_qubit_flip_locs = [x for x, y in data_qubit_flips.items() if y == 1]
        return data_qubit_flip_locs

    def _net_qubit_flips(
        self, total_flips: Dict[str, List[TQubit]]
    ) -> Dict[TQubitLoc, np.ndarray]:
        """Return all data qubits physical locations (no time dim) with a net flip

        Args:
            total_flips ({str:[node,]}): 
                key: flip type, either "X"/"Z"
                value: list of physical qubits (t,i,j) with a flip

        Returns:
            {(i,j),np.ndarray}:
                key: physical grid locaiton (i,j) of affected qubit
                value: net flip matrix on qubit
        """
        flips_x = total_flips["Z"]  # recall that the Z syndrome graph detects X Flips
        flips_z = total_flips["X"]  # recall that the X syndrome graph detects Z Flips
        flipsx = {flip: "X" for flip in flips_x if flip not in flips_z}
        flipsz = {flip: "Z" for flip in flips_z if flip not in flips_x}
        flipsy = {flip: "Y" for flip in flips_x if flip in flips_z}
        flips = {**flipsx, **flipsy, **flipsz}

        individual_flips: Dict[TQubitLoc, Dict[float, str]] = defaultdict(dict)

        for flip, error_key in flips.items():
            individual_flips[flip[1:]][flip[0]] = error_key

        physical_qubit_flips = {}
        for qubit_loc, flip_record in individual_flips.items():
            net_error = self.paulis["I"]
            # print("Physical Qubit: " + str(qubit_loc))
            for time, error in sorted(flip_record.items(), key=lambda item: item[0]):
                # print("Error: " + error + " at time: " + str(time))
                net_error = net_error.dot(self.paulis[error])
            physical_qubit_flips[qubit_loc] = net_error

        physical_qubit_flips = {
            x: y
            for x, y in physical_qubit_flips.items()
            if not np.array_equal(y, self.paulis["I"])
        }
        return physical_qubit_flips

    def _corrections(
        self, syndromes: Dict[str, List[TQubit]]
    ) -> Dict[TQubitLoc, np.ndarray]:
        """
        Args:
            syndromes ({str,[node,]}):
                key: syndrome_graph_key, either "X", "Z"
                value: activated syndrome nodes (t,i,j)

                Dictionary with syndromes["X"] containing tuples of the form 
                (t,i,j) where t specifies time and (i,j) specify position of 
                the X syndrome node changed from its value at (t-1,i,j), 
                and similarly, syndromes["Z"] for Z syndrome nodes. 

        Returns:
            net_flips ({(i,j):np.ndarray}): 
                dictionary with key representing physical (data) qubit
                and value representing the net error matrix on that data qubit. 
                e.g. key: (0,0), value: [[0,1],[1,0]] (X error)

        Additional Information:
            This method can be used to correct readout errors as shown in self.correct_readout. 
        """
        flips: Dict[str, List[TQubit]] = {}
        for syndrome_graph_key, syndrome_set in syndromes.items():
            if not syndrome_set:
                flips[syndrome_graph_key] = []
                continue
            error_graph, paths = self._make_error_graph(
                syndrome_set, syndrome_graph_key
            )  # TODO add option to use degeneracy weighting by setting err_prob
            matching_graph = self._make_matching_graph(error_graph, syndrome_graph_key)
            matches = self._run_mwpm(matching_graph, syndrome_graph_key)
            flips[syndrome_graph_key] = self._calculate_qubit_flips(
                matches, paths, syndrome_graph_key
            )
        net_flips = self._net_qubit_flips(flips)
        return net_flips

    def correct_readout(
        self,
        syndromes: Union[str, Dict[str, List[TQubit]]],
        logical_qubit_value: Optional[int] = None,
    ) -> int:
        """
        Args:
            readout: string like "1 00000000 00000000" representing "R S2 S1" (d=3, T=2) where 
            S1 is the first set of changed syndrome nodes (XOR'd with quiescent state syndrome measurements)
            S1 has the form: X3X2X1X0Z3Z2Z1Z0 in the case of d = 3. R represents the logical Z readout result.
        Returns:
            The most probable encoded value of the logical qubit.
        Additional Information:
            This method can be used to benchmark logical error rates, as well as perform fault tolerant readout.
        """
        if type(syndromes) == str:
            logical_qubit_value, syndromes = self._convert_string_to_nodes(
                str(syndromes)
            )
        syndromes = cast(Dict[str, List[TQubit]], syndromes)
        logical_qubit_value = cast(int, logical_qubit_value)
        # TODO is there a neater way to satisfy the type linter?

        # Logical Z readout will be performed with data qubits in the top row, this can be generalized later TODO
        net_flips = self._corrections(syndromes)
        top_row = [(0, i) for i in range(self.code_params["d"])]
        for top_data_qubit in top_row:
            if top_data_qubit in net_flips:
                flip_matrix = net_flips[top_data_qubit]
                error = self._find_pauli_matrix(flip_matrix)
                if error == "X" or error == "Y":
                    logical_qubit_value = (logical_qubit_value + 1) % 2
        return logical_qubit_value

    def _find_pauli_matrix(self, matrix: np.ndarray) -> str:
        for key, pauli in self.paulis.items():
            prod = matrix @ pauli
            if prod[0, 0] == 0:
                continue
            prod = prod / prod[0, 0]
            if np.all(prod == np.eye(2)):
                return key
        raise Exception("Not a Pauli Matrix")

    def _convert_string_to_nodes(
        self, readout_string: str
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
        chunks = readout_string.split(" ")

        int_syndromes = [int(x, base=2) for x in chunks[-1:0:-1]]
        xor_syndromes = [a ^ b for (a, b) in zip(int_syndromes, int_syndromes[1:])]

        X_syndromes = [(x & 0xF0) >> 4 for x in xor_syndromes]
        Z_syndromes = [x & 0xF for x in xor_syndromes]

        X = []
        for T, syndrome in enumerate(X_syndromes):
            for loc in range(4):
                if syndrome & 1 << loc:
                    X.append((float(T), -0.5 + loc, 0.5 + loc % 2))

        Z = []
        for T, syndrome in enumerate(Z_syndromes):
            for loc in range(4):
                if syndrome & 1 << loc:
                    Z.append((float(T), 0.5 + loc // 2, 0.5 + loc % 2 * 2 - loc // 2))

        return (
            int(chunks[0]),
            {"X": X, "Z": Z,},
        )

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
            x: plt.cm.plasma((y["time"] + 1) / self.code_params["T"])
            for x, y in G.nodes(data=True)
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


class XXZZGraphDecoder(VisualizationMixin, MWPMDecodingMixin, XXZZGraphDecoderBase):
    pass

