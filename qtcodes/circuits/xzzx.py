"""
XZZX Surface Code Encoder Classes
"""
from typing import Tuple, Optional
from qiskit import ClassicalRegister
from qiskit.circuit import Qubit

from qtcodes.circuits.base import _Stabilizer
from qtcodes.circuits.rotated_surface import _RotatedLattice, RotatedQubit

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

    def reset_x(self):
        """
        Initialize/reset to a logical |x+> state.
        """
        self.circ.reset(self.qregisters["data"])
        self.circ.h(self.qregisters["data"][1::2])  # H|0> = |+>

    def reset_z(self):
        """
        Initialize/reset to a logical |z+> state.
        """
        self.circ.reset(self.qregisters["data"])
        self.circ.h(self.qregisters["data"][0::2])  # H|0> = |+>

    def x(self) -> None:
        """
        Logical X operator on the qubit.
        Defined as the left-most column.
        """

        # Taking left-most column
        for i in range(0, int(self.params["num_data"]), int(self.params["d"])):
            if i % 2 == 1:
                self.circ.x(self.qregisters["data"][i])
            else:
                self.circ.z(self.qregisters["data"][i])
        self.circ.barrier()

    def z(self) -> None:
        """
        Logical Z operator on the qubit.
        Defined as the top-most row.
        """

        # Taking top-most row
        for i in range(int(self.params["d"])):
            if i % 2 == 0:
                self.circ.x(self.qregisters["data"][i])
            else:
                self.circ.z(self.qregisters["data"][i])
        self.circ.barrier()

    def x_c_if(self, classical: ClassicalRegister, val: int) -> None:
        """
        Classically conditioned logical X operator on the topological qubit.
        Defined as the left-most column.
        """

        # Taking left-most column
        for i in range(0, int(self.params["num_data"]), int(self.params["d"])):
            if i % 2 == 1:
                self.circ.x(self.qregisters["data"][i]).c_if(classical, val)
            else:
                self.circ.z(self.qregisters["data"][i]).c_if(classical, val)
        self.circ.barrier()

    def z_c_if(self, classical: ClassicalRegister, val: int) -> None:
        """
        Classically conditioned logical Z operator on the topological qubit.
        Defined as the top-most row.
        """

        # Taking top-most row
        for i in range(int(self.params["d"])):
            if i % 2 == 0:
                self.circ.x(self.qregisters["data"][i]).c_if(classical, val)
            else:
                self.circ.z(self.qregisters["data"][i]).c_if(classical, val)
        self.circ.barrier()

    def cx(self, control: Optional[Qubit] = None, target: Optional[Qubit] = None):
        """
        Logical CX Gate

        Args:
            control (Optional[Qubit]): If provided, then this gate will implement
                a logical x gate on this tqubit conditioned on source

            target (Optional[Qubit]): If provided, then this gate will implement
                a logical x gate on target conditioned on this tqubit
        """
        if control:
            # Taking left-most column
            for i in range(0, int(self.params["num_data"]), int(self.params["d"])):
                if i % 2 == 1:
                    self.circ.cx(control, self.qregisters["data"][i])
                else:
                    self.circ.cz(control, self.qregisters["data"][i])
            self.circ.barrier()
        elif target:
            self._readout_z_into_ancilla()
            self.circ.cx(self.qregisters["ancilla"], target)

    def _readout_x_into_ancilla(self) -> None:
        """
        Convenience method to read-out the
        logical-X projection into an ancilla qubit.
        Uses the left-most column.
        """

        self.circ.reset(self.qregisters["ancilla"])

        # Taking left-most column
        data_qubit_indxs = list(
            range(0, int(self.params["num_data"]), int(self.params["d"]))
        )

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

    def readout_x(self, readout_creg: Optional[ClassicalRegister] = None) -> None:
        """
        Convenience method to read-out the logical-X projection.
        Defined as the left-most column.
        """
        if not readout_creg:
            self.params["num_readout"] += 1
            creg_name = self.name + "_readout_" + str(self.params["num_readout"])
            readout = ClassicalRegister(1, name=creg_name)

            self.circ.add_register(readout)

            self.cregisters[creg_name] = readout
            readout_creg = self.cregisters[creg_name]
        self._readout_x_into_ancilla()
        self.circ.measure(self.qregisters["ancilla"], readout_creg)
        self.circ.barrier()

    def _readout_z_into_ancilla(self) -> None:
        """
        Convenience method to read-out the
        logical-Z projection into an ancilla qubit.
        Uses the top-most row.
        """

        self.circ.reset(self.qregisters["ancilla"])

        # Taking top-most row
        data_qubit_indxs = list(range(int(self.params["d"])))

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

    def readout_z(self, readout_creg: Optional[ClassicalRegister] = None) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        Defined as the top-most row.
        """
        if not readout_creg:
            self.params["num_readout"] += 1
            creg_name = self.name + "_readout_" + str(self.params["num_readout"])
            readout = ClassicalRegister(1, name=creg_name)

            self.circ.add_register(readout)

            self.cregisters[creg_name] = readout
            readout_creg = self.cregisters[creg_name]
        self._readout_z_into_ancilla()
        self.circ.measure(self.qregisters["ancilla"], readout_creg)
        self.circ.barrier()

    def lattice_readout_x(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of A stabilizer measurments,
        as well as a logical X readout.
        """

        self.params["num_lattice_readout"] += 1

        creg_name = (
            self.name + "_lattice_readout_" + str(self.params["num_lattice_readout"])
        )
        readout = ClassicalRegister(self.params["num_data"], name=creg_name)

        self.circ.add_register(readout)
        self.cregisters[creg_name] = readout

        # H|+> = |0>, H|-> = |1>
        # add H to measure along X for odd data qubits
        self.circ.h(self.qregisters["data"][1::2])

        self.circ.measure(self.qregisters["data"], self.cregisters[creg_name])
        self.circ.barrier()

    def lattice_readout_z(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of B stabilizer measurments,
        as well as a logical Z readout.
        """

        self.params["num_lattice_readout"] += 1
        creg_name = (
            self.name + "_lattice_readout_" + str(self.params["num_lattice_readout"])
        )
        readout = ClassicalRegister(self.params["num_data"], name=creg_name)

        self.circ.add_register(readout)
        self.cregisters[creg_name] = readout

        # H|+> = |0>, H|-> = |1>
        # add H to measure along X for even data qubits
        self.circ.h(self.qregisters["data"][0::2])

        self.circ.measure(self.qregisters["data"], self.cregisters[creg_name])
        self.circ.barrier()


class XZZXQubit(RotatedQubit):
    """
    A single logical surface code qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend TopologicalQubit which extends QuantumCircuit.
    """

    lattice_type = _XZZXLattice
