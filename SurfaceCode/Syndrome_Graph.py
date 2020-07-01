# -*- coding: utf-8 -*-
"""
This qiskit code was created on 29-JUN-20 1:18PM at IBM Hackatthon 2020 (Summer Jam)

@author: Shraddha Singh
"""

"""Generates syndrome graph for class GraphDecoder. This code has been modified 
using the code that is licensed under the Apache License, Version 2.0. You may
 obtain a copy of this license in the LICENSE.txt file in the root directory
of this source tree or at http://www.apache.org/licenses/LICENSE-2.0."""

import qiskit
from qiskit import QuantumRegister, ClassicalRegister
import copy
import warnings
import networkx as nx
import numpy as np
from qiskit import QuantumCircuit, execute
from random import randrange
from circuits import *
import matplotlib
try:
    from qiskit import Aer
    HAS_AER = True
except ImportError:
    from qiskit import BasicAer
    HAS_AER = False

    
class Syndrome():
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction code, and then run suitable fitters. 
    It imports objects from class circuits in surface_codes
    """

    def __init__(self, code, S=None):
        """
        Args:
            code (SurfaceCode(d,T)): The QEC Code object for which this decoder
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


    def _make_syndrome_graph(self):#syndrome graph X or Z
        """
        This method injects all possible Pauli errors into the circuit for
        ``code``.

        This is done by examining the qubits used in each gate of the
        circuit for a stored logical 0. A graph is then created with a node
        for each non-trivial syndrome element, and an edge between all such
        elements that can be created by the same error.
        """

        S = nx.Graph()

        qc = self.code.circuit['0']

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
                for error in ['x', 'z']:
                    temp_qc = copy.deepcopy(blank_qc)
                    temp_qc.name = str((j, qubit, error))
                    temp_qc.data = qc.data[0:j]
                    getattr(temp_qc, error)(qubit)
                    temp_qc.data += qc.data[j:depth + 1]
                    circuit_name[(j, qubit, error)] = temp_qc.name
                    error_circuit[temp_qc.name] = temp_qc

        if HAS_AER:
            simulator = Aer.get_backend('qasm_simulator')
        else:
            simulator = BasicAer.get_backend('qasm_simulator')

        job = execute(list(error_circuit.values()), simulator)

        for j in range(depth):
            qubits = qc.data[j][1]
            for qubit in qubits:
                for error in ['x', 'z']:

                    raw_results = {}
                    raw_results['0'] = job.result().get_counts(str((j, qubit, error)))
                    results = self.code.process_results(raw_results['0'])
                    
                    
                    
                    nodesX,nodesZ = self.code.extract_nodes(results)
                    for nodes in (nodesX,nodesZ):
                        for node in nodes:
                            print(node)
                            S.add_node(node)
                            for source in nodes:
                                for target in nodes:
                                    if source != target:
                                        S.add_edge(source, target, distance=1)

        return S
