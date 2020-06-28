# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 18:53:39 2020

@author: Shantanu Jha
"""

import networkx as nx
import numpy as np

from qiskit import QuantumCircuit, execute

class GraphDecoder():
	"""
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction code, and then run suitable decoders.
    """

    def __init__(self, code, S=None):
        """
        Args:
            code (RepitionCode): The QEC Code object for which this decoder
                will be used.
            S (networkx.Graph): Graph describing connectivity between syndrome
                elements. Will be generated automatically if not supplied.

        Additional information:
            The decoder for the supplied ``code`` is initialized by running
            ``_make_syndrome_graph()``. Since this process can take some
            time, it is also possible to load in a premade ``S``. However,
            if this was created for a differently defined ``code``, it won't
            work properly.
        """

        self.code = code

        if S:
            self.S = S
        else:
            self.S = self._make_syndrome_graph()

    def matching(self, error_graph):
    	matches = nx.max_weight_matching(error_graph, maxcardinality=True)
    	num_virtual = 0
    	for (source, target) in matches:
            num_virtual += E[source][target]['virtual']
        return matches, num_virtual
