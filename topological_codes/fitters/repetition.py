# -*- coding: utf-8 -*-
"""
Graph decoder for rep code
"""
from typing import Tuple, List, Dict

import retworkx as rx
from topological_codes.fitters.lattice_decoder import (
    LatticeGraphDecoder,
    TQubit,
)
from topological_codes.circuits.repetition import RepetitionQubit


class RepetitionGraphDecoder(LatticeGraphDecoder):
    encoder_type = RepetitionQubit

    def __init__(self, params: Dict) -> None:
        self.S["Z"] = rx.PyGraph(multigraph=False)
        self.node_map["Z"] = {}
        super().__init__(params)

    def _make_syndrome_graph(self) -> None:
        """
        Populates self.S["Z"] syndrome rx.PyGraph's with nodes specified by time and position.
        Args:
        Returns:
        """
        for vnode in self.virtual["Z"]:
            self.node_map["Z"][vnode] = self.S["Z"].add_node(vnode)

        edge_weight = 1
        for t in range(0, self.params["T"]):
            # real nodes
            for j in range(0, self.params["d"] - 1):
                node = (t, 0, j + 0.5)
                self.node_map["Z"][node] = self.S["Z"].add_node(node)

            # edges (real-real)
            for j in range(1, self.params["d"] - 1):
                left = (t, 0, j - 0.5)
                right = (t, 0, j + 0.5)
                self.S["Z"].add_edge(
                    self.node_map["Z"][left], self.node_map["Z"][right], edge_weight
                )

            # edges (real-virtual)
            self.S["Z"].add_edge(
                self.node_map["Z"][self.virtual["Z"][0]],
                self.node_map["Z"][(t, 0, 0.5)],
                edge_weight,
            )  # left
            self.S["Z"].add_edge(
                self.node_map["Z"][self.virtual["Z"][1]],
                self.node_map["Z"][(t, 0, self.params["d"] - 1.5)],
                edge_weight,
            )  # right

        # connect physical qubits in same location across subgraphs of adjacent times
        syndrome_nodes_t0 = [(t, x, y) for t, x, y in self.S["Z"].nodes() if t == 0]
        for node in syndrome_nodes_t0:
            space_label = (node[1], node[2])
            for t in range(0, self.params["T"] - 1):
                self.S["Z"].add_edge(
                    self.node_map["Z"][(t,) + space_label],
                    self.node_map["Z"][(t + 1,) + space_label],
                    edge_weight,
                )

    def _specify_virtual(self) -> Dict[str, List[TQubit]]:
        """
        Define coordinates of P virtual syndrome nodes (i.e. parity measurements to which we don't have access),
        Args:
        Returns:
            virtual (dictionary): where virtual["Z"] holds a list of tuples specifying virtual P syndrome nodes
        """
        virtual: Dict[str, List[TQubit]] = {}
        virtual["Z"] = []
        virtual["Z"].append((-1, 0, -0.5))
        virtual["Z"].append((-1, 0, self.params["d"] - 0.5))
        return virtual

    def _is_crossing_readout_path(
        self, match: Tuple[TQubit, TQubit], logical_readout_type: str
    ) -> bool:
        """
        Helper method that detects whether the match is crossing the readout path.

        Args:
            match (Tuple[TQubit, TQubit]): match in MWPM between two nodes
            logical_readout_type (str): e.g. "X"

        Returns:
            (bool): whether or not the match is crosing the readout path
        """
        raise NotImplementedError("Need to implement this...")
        return False
