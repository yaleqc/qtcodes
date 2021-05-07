"""
XZZX Surface Code Encoder Classes
"""
from typing import Tuple
from qiskit import ClassicalRegister
from .base import _Stabilizer
from .rotated_surface import _RotatedLattice, RotatedQubit

TQubit = Tuple[float, float, float]


class _XZZX(_Stabilizer):
    """
    XZZX-syndrome measurement for the rotated XZZX surface code.
    """

    def entangle(self) -> None:
        """
        Parity measurement on nearby data qubits.
        Going in a "Z" traversal direction.
        """
        syndrome = self.qubit_indices[0]
        top_l = self.qubit_indices[1]
        top_r = self.qubit_indices[2]
        bot_l = self.qubit_indices[3]
        bot_r = self.qubit_indices[4]

        if top_l:  # X
            self.circ.h(syndrome)
            self.circ.cx(syndrome, top_l)
            self.circ.h(syndrome)

        if top_r:  # Z
            self.circ.cx(top_r, syndrome)

        if bot_l:  # Z
            self.circ.cx(bot_l, syndrome)

        if bot_r:  # X
            self.circ.h(syndrome)
            self.circ.cx(syndrome, bot_r)
            self.circ.h(syndrome)


class _XZZXLattice(_RotatedLattice):
    """
    This class contains all the lattice geometry
    specifications regarding the XZZX Rotated Surface Code.
    """

    stabilizer_shortnames = {"mx": _XZZX, "mz": _XZZX}

    def logical_x_plus_reset(self):
        """
        Initialize/reset to a logical |x+> state.
        """
        self.circ.reset(self.qregisters["data"])
        self.circ.h(self.qregisters["data"][1::2])  # H|0> = |+>

    def logical_z_plus_reset(self):
        """
        Initialize/reset to a logical |z+> state.
        """
        self.circ.reset(self.qregisters["data"])
        self.circ.h(self.qregisters["data"][0::2])  # H|0> = |+>

    def logical_x(self) -> None:
        """
        Logical X operator on the qubit.
        Defined as the left-most column.
        """

        # Taking left-most column
        for i in range(0, self.params["num_data"], self.params["d"]):
            if i % 2 == 1:
                self.circ.x(self.qregisters["data"][i])
            else:
                self.circ.z(self.qregisters["data"][i])
        self.circ.barrier()

    def logical_z(self) -> None:
        """
        Logical Z operator on the qubit.
        Defined as the top-most row.
        """

        # Taking top-most row
        for i in range(self.params["d"]):
            if i % 2 == 0:
                self.circ.x(self.qregisters["data"][i])
            else:
                self.circ.z(self.qregisters["data"][i])
        self.circ.barrier()

    def readout_x(self) -> None:
        """
        Convenience method to read-out the logical-X projection.
        Defined as the left-most column.
        """
        self.params["num_readout"] += 1
        readout = ClassicalRegister(
            1, name=self.name + "_readout_" + str(self.params["num_readout"])
        )

        self.circ.add_register(readout)

        self.cregisters["readout"] = readout

        self.circ.reset(self.qregisters["ancilla"])

        # Taking left-most column
        data_qubit_indxs = list(range(0, self.params["num_data"], self.params["d"]))

        # X Readout
        x_readout_indxs = [i for i in data_qubit_indxs if i % 2 == 1]
        self.circ.h(self.qregisters["ancilla"])
        self.circ.cx(
            self.qregisters["ancilla"], self.qregisters["data"][x_readout_indxs]
        )
        self.circ.h(self.qregisters["ancilla"])

        # Z Readout
        z_readout_indxs = [i for i in data_qubit_indxs if i % 2 == 0]
        self.circ.cx(
            self.qregisters["data"][z_readout_indxs], self.qregisters["ancilla"]
        )
        self.circ.measure(self.qregisters["ancilla"], self.cregisters["readout"])
        self.circ.barrier()

    def readout_z(self) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        Defined as the top-most row.
        """

        self.params["num_readout"] += 1
        readout = ClassicalRegister(
            1, name=self.name + "_readout_" + str(self.params["num_readout"])
        )

        self.circ.add_register(readout)

        self.cregisters["readout"] = readout

        self.circ.reset(self.qregisters["ancilla"])

        # Taking top-most row
        data_qubit_indxs = list(range(self.params["d"]))

        # X Readout
        x_readout_indxs = [i for i in data_qubit_indxs if i % 2 == 0]
        self.circ.h(self.qregisters["ancilla"])
        self.circ.cx(
            self.qregisters["ancilla"], self.qregisters["data"][x_readout_indxs]
        )
        self.circ.h(self.qregisters["ancilla"])

        # Z Readout
        z_readout_indxs = [i for i in data_qubit_indxs if i % 2 == 1]
        self.circ.cx(
            self.qregisters["data"][z_readout_indxs], self.qregisters["ancilla"]
        )
        self.circ.measure(self.qregisters["ancilla"], self.cregisters["readout"])
        self.circ.barrier()

    def lattice_readout_x(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of A stabilizer measurments,
        as well as a logical X readout.
        """

        self.params["num_lattice_readout"] += 1
        readout = ClassicalRegister(
            self.params["num_data"],
            name=self.name
            + "_lattice_readout_"
            + str(self.params["num_lattice_readout"]),
        )

        self.circ.add_register(readout)
        self.cregisters["lattice_readout"] = readout

        # H|+> = |0>, H|-> = |1>
        # add H to measure along X for odd data qubits
        self.circ.h(self.qregisters["data"][1::2])

        self.circ.measure(self.qregisters["data"], self.cregisters["lattice_readout"])
        self.circ.barrier()

    def lattice_readout_z(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of B stabilizer measurments,
        as well as a logical Z readout.
        """

        self.params["num_lattice_readout"] += 1
        readout = ClassicalRegister(
            self.params["num_data"],
            name=self.name
            + "_lattice_readout_"
            + str(self.params["num_lattice_readout"]),
        )

        self.circ.add_register(readout)
        self.cregisters["lattice_readout"] = readout

        # H|+> = |0>, H|-> = |1>
        # add H to measure along X for even data qubits
        self.circ.h(self.qregisters["data"][0::2])

        self.circ.measure(self.qregisters["data"], self.cregisters["lattice_readout"])
        self.circ.barrier()


class XZZXQubit(RotatedQubit):
    """
    A single logical surface code qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend TopologicalQubit which extends QuantumCircuit.
    """

    lattice_type = _XZZXLattice
