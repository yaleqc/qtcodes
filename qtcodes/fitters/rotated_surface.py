# -*- coding: utf-8 -*-
"""
Graph decoder for surface codes
"""
from typing import Tuple, List, Dict

from qtcodes.circuits.xxzz import XXZZQubit
from qtcodes.fitters.lattice_decoder import (
    LatticeDecoder,
    TQubit,
    TQubitLoc,
)


class RotatedDecoder(LatticeDecoder):
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction surface code, and then run suitable decoders.
    """

    # TODO encoder is currently only used for string2node,
    # so we can use XXZZQubit, need to generalize
    encoder_type = XXZZQubit
    syndrome_graph_keys = ["X", "Z"]

    def _make_syndrome_graph(self) -> None:
        """
        Populates self.S["X"] and self.S["Z"] syndrome rx.PyGraph's
        with nodes specified by time and position.

        Args:
        Returns:
        """
        start_nodes = {"Z": (0.5, 0.5), "X": (0.5, 1.5)}
        for syndrome_graph_key in ["X", "Z"]:
            # subgraphs for each time step
            for t in range(0, self.params["T"]):
                start_node = start_nodes[syndrome_graph_key]
                node_label = (t,) + start_node
                self.node_map[syndrome_graph_key][node_label] = self.S[
                    syndrome_graph_key
                ].add_node(node_label)
                self._populate_syndrome_graph(
                    (t,) + start_node, t, [], syndrome_graph_key, 1
                )

            # connect physical qubits in same location across subgraphs of adjacent times
            syndrome_nodes_t0 = [
                (t, x, y) for t, x, y in self.S[syndrome_graph_key].nodes() if t == 0
            ]
            for node in syndrome_nodes_t0:
                space_label = (node[1], node[2])
                for t in range(0, self.params["T"] - 1):
                    self.S[syndrome_graph_key].add_edge(
                        self.node_map[syndrome_graph_key][(t,) + space_label],
                        self.node_map[syndrome_graph_key][(t + 1,) + space_label],
                        1,
                    )

    def _populate_syndrome_graph(
        self,
        current_node: TQubit,
        t: int,
        visited_nodes: List[TQubit],
        syndrome_graph_key: str,
        edge_weight: int = 1,
    ) -> None:
        """Recursive function to populate syndrome subgraph at time t with syndrome_graph_key X/Z.
        The current_node is connected to neighboring nodes without revisiting a node.

        Args:
            current_node ((t, x, y)): Current syndrome node to be connected with neighboring nodes.
            t (int): Current time, needed if current_node is a virtual node of the form (-1,i,j)
            visited_nodes ([(t, x, y),]): List of syndrome nodes which have already been traver.
            syndrome_graph_key (char): Which X/Z syndrome subgraph these nodes are from.
            edge_weight (float, optional): Weight of edge between two adjacent syndrome nodes.
                                           Defaults to 1.

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
            if target_node not in self.S[syndrome_graph_key].nodes():
                # add target_node to syndrome subgraph if it doesn't already exist
                self.node_map[syndrome_graph_key][target_node] = self.S[
                    syndrome_graph_key
                ].add_node(target_node)

            self.S[syndrome_graph_key].add_edge(
                self.node_map[syndrome_graph_key][current_node],
                self.node_map[syndrome_graph_key][target_node],
                edge_weight,
            )  # add edge between current_node and target_node

        # add virtual neighbors
        for target in virtual_neighbors:
            target_node = (
                -1,
            ) + target  # virtual target_node has time -1 with x and y coordinates from target
            if target_node not in self.S[syndrome_graph_key].nodes():
                # add virtual target_node to syndrome subgraph with
                # z coordinate (T-1)/2 for nice plotting,
                # if it doesn't already exist
                self.node_map[syndrome_graph_key][target_node] = self.S[
                    syndrome_graph_key
                ].add_node(target_node)
            self.S[syndrome_graph_key].add_edge(
                self.node_map[syndrome_graph_key][current_node],
                self.node_map[syndrome_graph_key][target_node],
                edge_weight,
            )  # add edge between current_node and virtual target_node

        # recursively traverse normal neighbors
        for target in normal_neighbors:
            self._populate_syndrome_graph(
                (t,) + target, t, visited_nodes, syndrome_graph_key, 1
            )

        # recursively traverse virtual neighbors
        for target in virtual_neighbors:
            self._populate_syndrome_graph(
                (-1,) + target, t, visited_nodes, syndrome_graph_key, 1
            )

    def _valid_syndrome(self, node: TQubitLoc, syndrome_graph_key: str) -> bool:
        """
        Checks whether a node is a syndrome node under our syndrome_graph_key,
        which is either X or Z.

        Args:
            node ((t, x, y)): Node in graph.
            syndrome_graph_key (char): Which X/Z syndrome subgraph these nodes are from.

        Returns:
            (bool): whether node is a syndrome node
        """
        i = node[0]
        j = node[1]
        if syndrome_graph_key == "Z":
            if i > 0 and i < self.params["d"] - 1 and j < self.params["d"] and j > -1:
                return True
            else:
                return False
        elif syndrome_graph_key == "X":
            if j > 0 and j < self.params["d"] - 1 and i < self.params["d"] and i > -1:
                return True
            else:
                return False
        else:
            raise ValueError("Please enter a valid syndrome_graph_key: X or Z")

    def _specify_virtual(self) -> Dict[str, List[TQubit]]:
        """
        Define coordinates of Z and X virtual syndrome nodes. Our convention is that Z
        virtual syndrome nodes are top/bottom and X virtual nodes are left/right.
        Args:
        Returns:
            virtual (dictionary): where virtual["X"] holds a list of tuples
            specifying virtual X syndrome nodes and equivalently for virtual["Z"]
        """
        virtual: Dict[str, List[TQubit]] = {}
        virtual["X"] = []
        virtual["Z"] = []
        for j in range(0, self.params["d"], 2):
            # Z virtual nodes
            virtual["Z"].append((-1, -0.5, j - 0.5))  # top
            virtual["Z"].append((-1, self.params["d"] - 0.5, j + 0.5))  # bottom

            # X virtual nodes
            virtual["X"].append((-1, j + 0.5, -0.5))  # left
            virtual["X"].append((-1, j - 0.5, self.params["d"] - 0.5))  # right
        return virtual

    def _is_crossing_readout_path(
        self, match: Tuple[TQubit, TQubit], logical_readout_type: str
    ) -> bool:
        """
        Helper method that detects whether the match is crossing the readout path.

        Args:
            match (Tuple[TQubit, TQubit]): match in MWPM between two nodes
            logical_readout_type (str): "X" or "Z"

        Returns:
            (bool): whether or not the match is crosing the readout path
        """

        # TODO Logical Z readout will be performed with data qubits in the top row,
        # this can be generalized later

        source, target = match
        if logical_readout_type == "Z":
            return (source[0] == -1 and source[1] == -0.5) or (
                target[0] == -1 and target[1] == -0.5
            )  # top
        elif logical_readout_type == "X":
            return (source[0] == -1 and source[2] == -0.5) or (
                target[0] == -1 and target[2] == -0.5
            )  # left
        else:
            raise ValueError("Please enter a valid logical_readout_type (X/Z).")
