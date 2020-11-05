# -*- coding: utf-8 -*-
"""
Generates circuits for quantum error correction with surface code patches.
"""
import copy
import warnings
from abc import ABC, abstractmethod

import numpy as np
import networkx as nx
from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, execute

try:
    from qiskit import Aer

    HAS_AER = True
except ImportError:
    from qiskit import BasicAer

    HAS_AER = False


class LatticeError(Exception):
    pass


class _Measure:
    def __init__(self, syndrome, top_l, top_r, bot_l, bot_r):
        self.syndrome = syndrome
        self.top_l = top_l
        self.top_r = top_r
        self.bot_l = bot_l
        self.bot_r = bot_r

    @abstractmethod
    def entangle(self, circ):
        pass


class _MeasureX(_Measure):
    def entangle(self, circ):
        """
        Traverse in reverse "Z" pattern
        """
        if (self.top_r and not self.top_l) or (self.bot_r and not self.bot_l):
            raise LatticeError("Inconsistent X syndrome connections")

        circ.h(self.syndrome)
        if self.top_r:
            circ.cx(self.syndrome, self.top_r)
            circ.cx(self.syndrome, self.top_l)
        if self.bot_r:
            circ.cx(self.syndrome, self.bot_r)
            circ.cx(self.syndrome, self.bot_l)
        circ.h(self.syndrome)


class _MeasureZ(_Measure):
    def entangle(self, circ):
        """
        Traverse in reverse "N" pattern
        """
        if (self.top_r and not self.bot_r) or (self.top_l and not self.bot_l):
            raise LatticeError("Inconsistent Z syndrome connections")

        if self.top_r:
            circ.cx(self.top_r, self.syndrome)
            circ.cx(self.bot_r, self.syndrome)
        if self.top_l:
            circ.cx(self.top_l, self.syndrome)
            circ.cx(self.bot_l, self.syndrome)


class RotatedSurfaceCodeLattice:
    def __init__(self, d, data_register, mx_register, mz_register):
        self.d = d
        self.measure_x = []
        self.measure_z = []

        per_row_x = (d ** 2 - 1) // 2 // (d + 1)
        per_row_z = (d ** 2 - 1) // 2 // (d - 1)
        for mx in mx_register:
            idx = mx.index
            row = idx // per_row_x
            offset = idx % per_row_x
            start = (row - 1) * d
            row_parity = row % 2

            if row == 0:  # First row
                top_l, top_r = None, None
                bot_l = data_register[idx * 2]
                bot_r = data_register[idx * 2 + 1]
            elif row == d:  # Last row
                bot_l, bot_r = None, None
                top_l = data_register[idx * 2 + 1]
                top_r = data_register[idx * 2 + 2]
            else:
                top_l = data_register[start + (offset * 2) + row_parity]
                top_r = data_register[start + (offset * 2) + row_parity + 1]
                bot_l = data_register[start + d + (offset * 2) + row_parity]
                bot_r = data_register[start + d + (offset * 2) + row_parity + 1]
            self.measure_x.append(_MeasureX(mx, top_l, top_r, bot_l, bot_r))

        for mz in mz_register:
            idx = mz.index
            row = idx // per_row_z
            offset = idx % per_row_z
            start = row * d
            row_parity = row % 2

            top_l = data_register[start + (offset * 2) - row_parity]
            top_r = data_register[start + (offset * 2) - row_parity + 1]
            bot_l = data_register[start + d + (offset * 2) - row_parity]
            bot_r = data_register[start + d + (offset * 2) - row_parity + 1]

            # Overwrite edge column syndromes
            if row_parity == 0 and offset == per_row_z:  # Last column
                top_r, bot_r = None, None
            elif row_parity == 1 and offset == 0:  # First column
                top_l, bot_l = None, None

            self.measure_z.append(_MeasureZ(mz, top_l, top_r, bot_l, bot_r))

    def entangle(self, circ):
        for syndrome in (self.measure_x, self.measure_z):
            for ancilla in syndrome:
                ancilla.entangle(circ)

    def parse_readout(self, readout_string):
        syn_len = (d ** 2 - 1) // 2
        chunks = readout_string.split(" ")

        int_syndromes = [int(x, base=2) for x in chunks[-1:0:-1]]
        xor_syndromes = [a ^ b for (a, b) in zip(int_syndromes, int_syndromes[1:])]

        mask_Z = "1" * syn_len
        mask_X = mask_Z + "0" * syn_len
        X_syndromes = [(x & int(mask_X, base=2)) >> syn_len for x in xor_syndromes]
        Z_syndromes = [x & int(mask_Z, base=2) for x in xor_syndromes]

        X = []
        for T, syndrome in enumerate(X_syndromes):
            for loc in range(syn_len):
                if syndrome & 1 << loc:
                    X.append((T, -0.5 + loc, 0.5 + loc % 2))

        Z = []
        for T, syndrome in enumerate(Z_syndromes):
            for loc in range(syn_len):
                if syndrome & 1 << loc:
                    Z.append((T, 0.5 + loc // 2, 0.5 + loc % 2 * 2 - loc // 2))

        return (
            int(chunks[0]),
            {"X": X, "Z": Z,},
        )


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
        self.T = T
        self.num_data = d ** 2
        self.num_syn = (d ** 2 - 1) // 2

        # These registers can be used as pointers across all of the circuits we create
        self.data = QuantumRegister(self.num_data, "data")
        self.mx = QuantumRegister(self.num_syn, "mx")
        self.mz = QuantumRegister(self.num_syn, "mz")
        self.measurements = [
            ClassicalRegister(self.num_syn * 2, name="c{}".format(i + 1))
            for i in range(T + 1)  # First round is quiescent projection
        ]

        # Generate the lattice connections within the convention we are using
        self.lattice = RotatedSurfaceCodeLattice(d, self.data, self.mx, self.mz)

    @property
    def zero(self):
        circ = QuantumCircuit(self.data, self.mz, self.mx, *self.measurements)

        for i in range(self.T + 1):
            self.lattice.entangle(circ)
            circ.barrier()
            circ.measure(self.mz, self.measurements[i][0 : self.num_syn])
            circ.measure(self.mx, self.measurements[i][self.num_syn : self.num_syn * 2])
            circ.reset(self.mz)
            circ.reset(self.mx)
            circ.barrier()

        return circ
