# -*- coding: utf-8 -*-

"""Generates circuits for quantum error correction with surface code patches."""

import qiskit
from qiskit import QuantumRegister, ClassicalRegister
import copy
import warnings
import networkx as nx
import numpy as np
from qiskit import QuantumCircuit, execute
from random import randrange
import matplotlib

try:
    from qiskit import Aer

    HAS_AER = True
except ImportError:
    from qiskit import BasicAer

    HAS_AER = False


class SurfaceCode:
    """
    Implementation of a distance d surface code, implemented over
    T syndrome measurement rounds.
    """

    def __init__(self, d, T):
        """
        Creates the circuits corresponding to a logical 0 encoded
        using a surface code with X and Z stabilizers.
        Args:
            d (int): Number of physical "data" qubits. Only odd d's allowed
            T (int): Number of rounds of ancilla-assisted syndrome measurement. Normally T=d
        Additional information:
            No measurements are added to the circuit if `T=0`. Otherwise
            `T` rounds are added, followed by measurement of the code
            qubits (corresponding to a logical measurement and final
            syndrome measurement round)
            This circuit is for "rotated lattices" i.e. it requires
            d**2 data qubits and d**2-1 syndrome qubits only. Hence,
            d=odd allows equal number of Z and X stabilizer mesaurments.
        """
        self.d = d
        self.T = 0
        self.data = QuantumRegister(d ** 2, "data")
        self.ancilla = QuantumRegister((d ** 2 - 1), "ancilla")
        self.qubit_registers = {"data", "ancilla"}
        self.output = []
        self.circuit = {}
        self.c_output = ClassicalRegister(d ** 2, "c_output")

        """This code creates circuits only for log='0' but it can be easily
        modified to accomodate circuuit for log='1' """

        for log in ["0", "1"]:
            self.circuit[log] = QuantumCircuit(self.ancilla, self.data, name=log)

        #        self._preparation() to be included to create log='1'

        for _ in range(T - 1):
            self.syndrome_measurement()

        if T != 0:
            self.syndrome_measurement(reset=False)
            self.readout()

    def get_circuit_list(self):
        """
        Returns:
            circuit_list: self.circuit as a list, with
            circuit_list[0] = circuit['0']
            circuit_list[1] = circuit['1']
        """
        circuit_list = [self.circuit[log] for log in ["0", "1"]]
        return circuit_list

    """It assigns vertices to qubits on a 2D graph, where,
    data qubits are on the x lines and, syndrome qubits are
    on the 0.5+x lines, in the cartesian coordinate system, where x is an integer."""

    def lattice(self):
        d = self.d
        data_string = nx.Graph()
        syndrome_string = nx.Graph()
        for i in range(0, d):
            for j in range(0, d):
                data_string.add_node((i, j))
        for k in range(0, d, 1):
            for i in range(0, d + 1, 1):
                for j in range(0, d + 1, 1):
                    if (i + j) % 2 != 0:
                        if ((i % 2 == 0) and j != d) or ((i % 2 == 1) and (j != 0)):
                            syndrome_string.add_node(((2 * i - 1) / 2, (2 * j - 1) / 2))
                    else:
                        if ((j % 2 == 0) and i != 0) or ((j % 2 == 1) and (i != d)):
                            syndrome_string.add_node(((2 * i - 1) / 2, (2 * j - 1) / 2))

        syn_ind = list(syndrome_string.nodes)
        data_ind = list(data_string.nodes)
        return (syn_ind, data_ind)

    """List of nodes on the 2D lattice graph returned"""
    # def x """to be included to execute self._preparation for log='1' """

    # def _preparation(self):
    """ prapares log '1' from log '0' circuit by applying x logical to the lattice"""

    def connection(self):
        """
        Determines the order of syndrome measurements between data qubits and syndrome qubits.
        We follow the ZN rule here to avoid hook error as described by [https://doi.org/10.1063/1.1499754]
        where Z stabilisers are arranged in 'Z' pattern and X stabilizers in 'N' pattern.
        Refer to the diagram in readme to get the refrence.
        """
        syn_index, data_index = self.lattice()

        order = []
        for i in range(self.d ** 2 - 1):
            d = data_index
            r = syn_index[i][0]
            c = syn_index[i][1]

            def get_index(j):
                for i in range(len(data_index)):
                    if data_index[i] == j:
                        return i

            new = []
            new.append((r, c))
            if r == -0.5:  # top semicircile
                new.append(-1)
                new.append(get_index((r + 0.5, c - 0.5)))
                new.append(-1)
                new.append(get_index((r + 0.5, c + 0.5)))
            elif c == -0.5:  # left semicircle
                new.append(-1)
                new.append(get_index((r - 0.5, c + 0.5)))
                new.append(-1)
                new.append(get_index((r + 0.5, c + 0.5)))

            elif r == self.d - 0.5:  # bottom semicircle

                new.append(get_index((r - 0.5, c - 0.5)))
                new.append(-1)
                new.append(get_index((r - 0.5, c + 0.5)))
                new.append(-1)

            elif c == self.d - 0.5:  # right semicircle
                new.append(get_index((r - 0.5, c - 0.5)))
                new.append(-1)
                new.append(get_index((r + 0.5, c - 0.5)))
                new.append(-1)
            else:
                if (r + c) % 2 == 0:  # square patches
                    new.append(get_index((r - 0.5, c - 0.5)))
                    new.append(get_index((r + 0.5, c - 0.5)))
                    new.append(get_index((r - 0.5, c + 0.5)))
                    new.append(get_index((r + 0.5, c + 0.5)))
                else:
                    new.append(get_index((r - 0.5, c - 0.5)))
                    new.append(get_index((r - 0.5, c + 0.5)))
                    new.append(get_index((r + 0.5, c - 0.5)))
                    new.append(get_index((r + 0.5, c + 0.5)))
            order.append(new)
        return order

    def syndrome_measurement(self, reset=True, barrier=True):
        """
            Application of a syndrome measurement round.
            Args:
                reset (bool): If set to true add a boolean at the end of each round
                barrier (bool): Boolean denoting whether to include a barrier at the end.
                A barrier is included after every round of 'j' which passes through layers of
                cx to be done, because the order should not be disturbed else the stabilizers
                will not be executed since Z and X on the same qubit do not commute. Thus,
                we end up flipping the sign of some stabilizers.
            """
        self.output.append(
            ClassicalRegister((self.d ** 2 - 1), "round_" + str(self.T) + "ancilla")
        )

        for log in ["0", "1"]:
            self.circuit[log].add_register(self.output[-1])
            order = self.connection()
            for j in range(1, 5):
                for i in range(len(order)):
                    k = self.data[order[i][j]]
                    l = self.ancilla[i]
                    if (order[i][0][0] + order[i][0][1]) % 2 == 0:  # Xstabilizer
                        if j == 1:
                            self.circuit[log].h(l)
                        if order[i][j] != -1:
                            self.circuit[log].cx(l, k)
                        if j == 4:
                            self.circuit[log].h(l)
                    else:  # Xstabilizer
                        if order[i][j] != -1:
                            self.circuit[log].cx(k, l)
                if barrier:
                    self.circuit[log].barrier()

            for j in range(self.d ** 2 - 1):
                if (order[j][0][0] + order[j][0][1]) % 2 == 1:  # Z
                    self.circuit[log].measure(self.ancilla[j], self.output[self.T][j])
                if reset:
                    self.circuit[log].reset(self.ancilla[j])

            self.circuit[log].barrier()

            for j in range(self.d ** 2 - 1):
                if (order[j][0][0] + order[j][0][1]) % 2 == 0:  # X
                    self.circuit[log].measure(self.ancilla[j], self.output[self.T][j])
                if reset:
                    self.circuit[log].reset(self.ancilla[j])

        self.T += 1

    def readout(self):
        """
        Readout of all code qubits, which corresponds to a logical measurement
        as well as allowing for a measurement of the syndrome to be inferred.
        """
        for log in ["0", "1"]:
            self.circuit[log].add_register(self.c_output)
            for i in range(self.d ** 2):
                self.circuit[log].measure(self.data[i], self.c_output[i])

    def process_results(self, raw_results):
        """
        Args:
            raw_results (dict): A dictionary whose keys are logical values,
                and whose values are standard counts dictionaries, (as
                obtained from the `get_counts` method of a ``qiskit.Result``
                object).
        Returns:
            syn: d+1 dimensional array where 0th array stores qubit readouts
            while the subsequesnt rows store the results from measurement rounds
            as required for extraction of nodes with errors to be sent to the decoder
        Additional information:
            The circuits must be executed outside of this class, so that
            their is full freedom to compile, choose a backend, use a
            noise model, etc. The results from these executions should then
            be used to create the input for this method.
        """
        results = []
        results = list(max(raw_results, key=raw_results.get))

        syn = []
        new = []
        for i in results:
            for j in range(len(i)):
                if i[j] != " ":
                    new.append(int(i[j]))
                else:
                    syn.append(new)
                    new = []
        syn.append(new)

        return syn

    def extract_nodes(self, syn_meas_results):
        """Extracts node locations of qubits which flipped in
        consecutive rounds (stored as (k,i,j)) and the data qubits which were flipped
        during readout (stored as (-2,i,j)). Here k spans range(0,d-1,1)
        Z syndrome nodes and Z logical data qubit nodes (see figure) in error_nodesZ
        and we do the same for X stabilizers and X logical qubits in error_nodesX.
        Note that arrays are reversed in terms of syndrome rounds, when compared to
        syn_meas_results
        """
        processed_results = []
        new = []
        for j in syn_meas_results[0]:
            new.append(j)
        processed_results.append(new)
        new = []
        for j in syn_meas_results[len(syn_meas_results) - 1]:
            new.append(j)
        processed_results.append(new)

        for i in range(len(syn_meas_results) - 2, 0, -1):
            new = []
            for j in range(0, len(syn_meas_results[i])):
                new.append((syn_meas_results[i][j] + syn_meas_results[i + 1][j]) % 2)
            processed_results.append(new)

        syn, dat = self.lattice()
        error_nodesX = []
        error_nodesZ = []

        # first_row = processed_result[0][:self.d]
        # last_row = processed_result[0][-self.d - 1:-1]

        # left_col = processed_result[0][::self.d]
        # right_col = processed_result[0][self.d-1:-1:self.d]

        # if sum(first_row) % 2 == 1 or sum(last_row) % 2 == 1:
        #     for node in dat[:self.d]:
        #         # Append virtual node
        #         if node[1] == 0:
        #             error_nodesZ.append((-1, node[0] - 0.5, node[1] - 0.5))
        #         else:
        #             error_nodesZ.append((-1, node[0] - 0.5, node[1] + 0.5))

        #     for node in dat[-self.d - 1:-1]:
        #         if node[1] == self.d - 1:
        #             error_nodesZ.append((-1, node[0] + 0.5, node[1] + 0.5))
        #         else:
        #             error_nodesZ.append((-1, node[0] + 0.5, node[1] - 0.5))

        # if sum(left_col) % 2 == 1 or sum(right_col) % 2 == 1:
        #     for node in dat[::self.d]:
        #         error_nodesX.append((-2, node[0], node[1]))
        #     for node in dat[self.d-1:-1:self.d]:
        #         error_nodesX.append((-2, node[0], node[1]))

        for i in range(1, len(processed_results)):
            for j in range(len(processed_results[i])):

                if processed_results[i][j] == 1:

                    if (syn[j][0] + syn[j][1]) % 2 == 0:
                        error_nodesX.append((i - 1, syn[j][0], syn[j][1]))
                    else:
                        error_nodesZ.append((i - 1, syn[j][0], syn[j][1]))
        return error_nodesX, error_nodesZ
