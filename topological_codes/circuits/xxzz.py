"""
XXZZ Surface Code Encoder Classes
"""
from typing import Tuple
from qiskit import ClassicalRegister
from topological_codes.circuits.base import (
    _Stabilizer,
    LatticeError,
)
from topological_codes.circuits.rotated_surface import _RotatedLattice, RotatedQubit


TQubit = Tuple[float, float, float]


class _XXXX(_Stabilizer):
    """
    X-syndrome plaquette of the rotated CSS (XXXX/ZZZZ) surface code.
    """

    def entangle(self) -> None:
        """
        Traverse in reverse "Z" pattern
        """
        syndrome = self.qubit_indices[0]
        top_l = self.qubit_indices[1]
        top_r = self.qubit_indices[2]
        bot_l = self.qubit_indices[3]
        bot_r = self.qubit_indices[4]

        if (top_r and not top_l) or (bot_r and not bot_l):
            raise LatticeError("Inconsistent X syndrome connections")

        self.circ.h(syndrome)
        if top_r:
            self.circ.cx(syndrome, top_r)
            self.circ.cx(syndrome, top_l)
        if bot_r:
            self.circ.cx(syndrome, bot_r)
            self.circ.cx(syndrome, bot_l)
        self.circ.h(syndrome)


class _ZZZZ(_Stabilizer):
    """
    Z-syndrome plaquette of the rotated surface code.
    """

    def entangle(self) -> None:
        """
        Traverse in reverse "N" pattern
        """
        syndrome = self.qubit_indices[0]
        top_l = self.qubit_indices[1]
        top_r = self.qubit_indices[2]
        bot_l = self.qubit_indices[3]
        bot_r = self.qubit_indices[4]

        if (top_r and not bot_r) or (top_l and not bot_l):
            raise LatticeError("Inconsistent Z syndrome connections")

        if top_r:
            self.circ.cx(top_r, syndrome)
            self.circ.cx(bot_r, syndrome)
        if top_l:
            self.circ.cx(top_l, syndrome)
            self.circ.cx(bot_l, syndrome)


class _XXZZLattice(_RotatedLattice):
    """
    This class contains all the lattice geometry
    specifications regarding the XXZZ (CSS) Rotated Surface Code.
    """

    stabilizer_shortnames = {"mx": _XXXX, "mz": _ZZZZ}

    def logical_plus_x_reset(self) -> None:
        """
        Initialize/reset to a logical |x+> state.
        """
        self.circ.reset(self.qregisters["data"])
        self.circ.h(self.qregisters["data"])
        self.circ.barrier()

    def logical_plus_z_reset(self) -> None:
        """
        Initialize/reset to a logical |z+> state.
        """
        self.circ.reset(self.qregisters["data"])
        self.circ.barrier()

    def logical_x(self) -> None:
        """
        Logical X operator on the qubit.
        Uses the left-most column.
        """
        for i in range(0, self.params["num_data"], self.params["d"]):
            self.circ.x(self.qregisters["data"][i])
        self.circ.barrier()

    def logical_z(self) -> None:
        """
        Logical Z operator on the qubit.
        Uses the top-most row.
        """
        for i in range(self.params["d"]):
            self.circ.z(self.qregisters["data"][i])
        self.circ.barrier()

    def readout_x(self) -> None:
        """
        Convenience method to read-out the logical-X projection.
        Uses the left-most column.
        """
        self.params["num_readout"] += 1
        readout = ClassicalRegister(
            1, name=self.name + "_readout_" + str(self.params["num_readout"])
        )

        self.circ.add_register(readout)

        self.cregisters["readout"] = readout

        self.circ.reset(self.qregisters["ancilla"])
        self.circ.h(self.qregisters["ancilla"])
        for i in range(0, self.params["num_data"], self.params["d"]):
            self.circ.cx(self.qregisters["ancilla"], self.qregisters["data"][i])
        self.circ.h(self.qregisters["ancilla"])
        self.circ.measure(self.qregisters["ancilla"], self.cregisters["readout"])
        self.circ.barrier()

    def readout_z(self) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        Uses the top-most row.
        """
        self.params["num_readout"] += 1
        readout = ClassicalRegister(
            1, name=self.name + "_readout_" + str(self.params["num_readout"])
        )

        self.circ.add_register(readout)

        self.cregisters["readout"] = readout

        self.circ.reset(self.qregisters["ancilla"])
        for i in range(self.params["d"]):
            self.circ.cx(self.qregisters["data"][i], self.qregisters["ancilla"])
        self.circ.measure(self.qregisters["ancilla"], self.cregisters["readout"])
        self.circ.barrier()

    def lattice_readout_x(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of X stabilizer measurments,
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
        self.circ.h(self.qregisters["data"])
        self.circ.measure(self.qregisters["data"], self.cregisters["lattice_readout"])
        self.circ.barrier()

    def lattice_readout_z(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of Z stabilizer measurments,
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

        self.circ.measure(self.qregisters["data"], self.cregisters["lattice_readout"])
        self.circ.barrier()


class XXZZQubit(RotatedQubit):
    """
    A single logical surface code qubit.
    """

    lattice_type = _XXZZLattice
