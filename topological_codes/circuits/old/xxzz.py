# -*- coding: utf-8 -*-
"""
Generates circuits for quantum error correction with surface code patches.
"""
import copy
import warnings
from abc import abstractmethod

import numpy as np
import networkx as nx
from typing import TypeVar, Tuple, Dict, List
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, execute

from topological_codes import LatticeError, TopologicalQubit, TopologicalLattice

try:
    from qiskit import Aer

    HAS_AER = True
except ImportError:
    from qiskit import BasicAer

    HAS_AER = False

TQubit = Tuple[float, float, float]


class _Face:
    """
    Abstract class for a single *face* of the surface code, which is described
    by a syndrome qubit and the four data qubits that surround it. If the face
    exists on an edge of the surface code lattice, some of the data qubits may
    be None.
    """

    def __init__(self, syndrome, top_l, top_r, bot_l, bot_r):
        """
        Initialize a face by passing in the single QubitRegisters of the circuit
        in which this face will be embedded.
        """
        self.syndrome = syndrome
        self.top_l = top_l
        self.top_r = top_r
        self.bot_l = bot_l
        self.bot_r = bot_r

    @abstractmethod
    def entangle(self, circ):
        pass


class _FaceX(_Face):
    """
    X-syndrome face of the rotated surface code.
    """

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


class _FaceZ(_Face):
    """
    Z-syndrome face of the rotated surface code.
    """

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


class RotatedSurfaceCodeLattice(TopologicalLattice):
    """
    This class is essentially a helper class that translates between the lattice
    abstraction of the surface code and the physical qubits in the circuit. This
    promotes a clean separation of concerns and wraps the implementation details
    of the lattice + code entanglement ordering.
    """

    def __init__(
        self,
        d: int,
        data_register: QuantumRegister,
        mx_register: QuantumRegister,
        mz_register: QuantumRegister,
    ) -> None:
        """
        Initializes an instance of the rotated surface code lattice with our
        chosen layout and numbering.

        Args:
            d (int): surface code distance
            data_register (QuantumRegister): grouped register of all data qubits
            mx_register (QuantumRegister): grouped register of all measure-x qubits
            mz_register (QuantumRegister): grouped register of all measure-z qubits
        """
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
            self.measure_x.append(_FaceX(mx, top_l, top_r, bot_l, bot_r))

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
            if row_parity == 0 and offset == per_row_z - 1:  # Last column
                top_r, bot_r = None, None
            elif row_parity == 1 and offset == 0:  # First column
                top_l, bot_l = None, None

            self.measure_z.append(_FaceZ(mz, top_l, top_r, bot_l, bot_r))
        super().__init__(d, data_register, mx_register, mz_register)

    def entangle(self, circ: QuantumCircuit) -> None:
        """
        Entangles the entire surface code by calling the entangle method of each
        syndrome face. Within a face, order is determined by the delegated
        method. Order between faces should not matter here, since the projection
        will be determined by the measurement order.
        """
        syn_len = (self.d ** 2 - 1) // 2
        for loc in range(syn_len):
            self.measure_x[loc].entangle(circ)
            self.measure_z[loc].entangle(circ)
            circ.barrier()

    def parse_readout(self, readout_string: str) -> Tuple[int, Dict[str, List[TQubit]]]:
        """
        Helper method to turn a result string (e.g. 1 10100000 10010000) into an
        appropriate logical readout value and XOR-ed syndrome locations
        according to our grid coordinate convention.
        """
        syn_len = (self.d ** 2 - 1) // 2
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


class XXZZQubit(TopologicalQubit):
    """
    A single logical surface code qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend QuantumCircuit.
    """

    def __init__(self, d: int, *args, **kwargs) -> None:
        """
        Initializes a new QuantumCircuit for this logical qubit and calculates
        the underlying surface code lattice ordering.
        
        Args:
            d (int): Number of physical "data" qubits. Only odd d is possible!
        """
        if d % 2 != 1:
            raise ArgumentError("Surface code distance must be odd!")
        self.__d = d
        self.__T = 0
        self.__num_data = d ** 2
        self.__num_syn = (d ** 2 - 1) // 2

        self.__data = QuantumRegister(self.__num_data, "data")
        self.__mz = QuantumRegister(self.__num_syn, "mz")
        self.__mx = QuantumRegister(self.__num_syn, "mx")

        # We implement and assume the rotated lattice only, but this can be
        # imagined to accept other layouts in the future.
        self.__lattice = RotatedSurfaceCodeLattice(d, self.__data, self.__mx, self.__mz)

        # Spare ancilla (e.g. for readout)
        self.__ancilla = QuantumRegister(1, name="ancilla")
        super().__init__(self.__data, self.__mz, self.__mx, self.__ancilla)

    def stabilize(self) -> None:
        """
        Run a single round of stabilization (entangle and measure).
        """
        syndrome_readouts = ClassicalRegister(
            self.__num_syn * 2, name="c{}".format(self.__T)
        )
        self.add_register(syndrome_readouts)
        self.__T += 1

        self.__lattice.entangle(self)
        self.measure(self.__mz, syndrome_readouts[0 : self.__num_syn])
        self.measure(self.__mx, syndrome_readouts[self.__num_syn : self.__num_syn * 2])
        self.reset(self.__mz)
        self.reset(self.__mx)
        self.barrier()

    def identity(self) -> None:
        """
        Inserts an identity on the data and syndrome qubits. This is a hack to
        create an isolated error model.
        """
        [
            self.id(x)
            for register in (self.__data, self.__mz, self.__mx)
            for x in register
        ]
        self.barrier()

    def identity_data(self) -> None:
        """
        Inserts an identity on the data qubits only. This is a hack to create an
        isolated error model.
        """
        [self.id(x) for register in (self.__data,) for x in register]
        self.barrier()

    def hadamard_reset(self) -> None:
        """
        A hack to initialize a + and - logical qubit for now...
        """
        [self.reset(x) for x in self.__data]
        [self.h(x) for x in self.__data]
        self.barrier()

    def logical_x(self) -> None:
        """
        Logical X operator on the qubit.
        """
        for i in range(0, self.__num_data, self.__d):
            self.x(self.__data[i])
        self.barrier()

    def logical_z(self) -> None:
        """
        Logical Z operator on the qubit.
        """
        for i in range(self.__d):
            self.z(self.__data[i])
        self.barrier()

    def readout_z(self) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        """
        readout = ClassicalRegister(1, name="readout")
        self.add_register(readout)

        self.reset(self.__ancilla)
        for i in range(self.__d):
            self.cx(self.__data[i], self.__ancilla)
        self.measure(self.__ancilla, readout)
        self.barrier()

    def readout_x(self) -> None:
        """
        Convenience method to read-out the logical-X projection.
        """
        readout = ClassicalRegister(1, name="readout")
        self.add_register(readout)

        self.reset(self.__ancilla)
        self.h(self.__ancilla)
        for i in range(0, self.__num_data, self.__d):
            self.cx(self.__ancilla, self.__data[i])
        self.h(self.__ancilla)
        self.measure(self.__ancilla, readout)
        self.barrier()

    def parse_readout(self, readout_string: str) -> Tuple[int, Dict[str, List[TQubit]]]:
        return self.__lattice.parse_readout(readout_string)
