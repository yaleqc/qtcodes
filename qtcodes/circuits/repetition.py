"""
Repetition Code Encoder Classes
"""
from typing import Dict, List, Tuple, Optional, Type
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Qubit

from qtcodes.circuits.base import (
    _Stabilizer,
    _TopologicalLattice,
    TopologicalQubit,
)


TQubit = Tuple[float, float, float]


class _Parity(_Stabilizer):
    """
    Parity syndrome measurement for Repetition code.
    """

    def entangle(self) -> None:
        """
        Parity measurement on nearby data qubits.
        """
        syndrome = self.qubit_indices[0]
        left = self.qubit_indices[1]
        right = self.qubit_indices[2]

        self.circ.cx(left, syndrome)
        self.circ.cx(right, syndrome)


class _RepetitionLattice(_TopologicalLattice):
    """
    This class contains all the lattice geometry specifications regarding the Repetition Code.

    E.g.
    X
    |
    X
    |
    X
    |
    X
    |
    X

    d=5 Rep Code
    """

    def __init__(self, params: Dict[str, float], name: str, circ: QuantumCircuit):
        """
        Initializes this Topological Lattice class.

        Args:
            params (Dict[str,int]):
                Contains params such as d, where d is the number of
                physical "data" qubits lining a row or column of the lattice.
            name (str):
                Useful when combining multiple TopologicalQubits together.
                Prepended to all registers.
            circ (QuantumCircuit):
                QuantumCircuit on top of which the topological qubit is built.
                This is often shared amongst multiple TQubits.
        """
        self.geometry: Dict[str, List[List[int]]] = {}
        super().__init__(params, name, circ)

    def _params_validate_and_generate(self) -> None:
        """
        Validate and generate params.

        E.g.
        self.params["num_syn"] = params["d"] - 1
        """
        # default params
        if "d" not in self.params:
            self.params["d"] = 3

        # calculated params
        self.params["T"] = -1  # -1 until a stabilizer round is added!
        self.params["num_readout"] = -1  # -1 until a logical readout is performed!
        self.params[
            "num_lattice_readout"
        ] = -1  # -1 until a lattice readout is performed!
        self.params["num_data"] = self.params["d"]
        self.params["num_syn"] = self.params["d"] - 1

    def _gen_registers(self) -> None:
        """
        Implement this method to create quantum and classical registers.

        E.g.
        qregisters["data"] = QuantumRegister(params["num_data"], name=name + "_data")
        """
        self.qregisters["data"] = QuantumRegister(
            self.params["num_data"], name=self.name + "_data"
        )
        self.qregisters["mz"] = QuantumRegister(
            self.params["num_syn"], name=self.name + "_mp"
        )
        self.qregisters["ancilla"] = QuantumRegister(1, name=self.name + "_ancilla")

    def _set_geometry(self):
        """
        Construct the lattice geometry for reuse across this class.

        Returns:
            geometry (Dict[str, List[List[int]]]):
                key: syndrome/plaquette type
                value: List of lists of qubit indices comprising one plaquette.
        """
        geometry = {"mz": []}

        for i in range(int(self.params["num_syn"])):
            syn = i
            left = i
            right = i + 1
            geometry["mz"].append([syn, left, right])

        self.geometry = geometry

    def _gen_qubit_indices_and_stabilizers(
        self,
    ) -> Tuple[List[List[Qubit]], List[Type[_Parity]]]:
        """
        Generates lattice blueprint for rep code with our
        chosen layout and numbering.

        Returns:
            qubit_indices (List[List[Qubit]]):
                List of lists of Qubits that comprise each plaquette.

            stabilizers (List[_Stabilizer]):
                List of stabilizers for each plaquette.
        """
        self._set_geometry()

        qubit_indices = []
        stabilizers = []
        for stabilizer, idx_lists in self.geometry.items():
            stabilizer_cls = _Parity
            for idx_list in idx_lists:
                syn = self.qregisters[stabilizer][idx_list[0]]
                plaquette = [
                    self.qregisters["data"][idx] if idx is not None else None
                    for idx in idx_list[1:]
                ]
                plaquette = [syn,] + plaquette
                qubit_indices.append(plaquette)
                stabilizers.append(stabilizer_cls)
        return qubit_indices, stabilizers

    def extract_final_stabilizer_and_logical_readout_z(
        self, final_readout_string: str
    ) -> Tuple[int, str]:
        """
        Extract final Parity syndrome measurements and logical Z readout from
        lattice readout along the Parity syndrome graph.

        Args:
            final_readout_string (str):
                readout string of length equal to the number of data qubits
                contains the readout values of each data qubit measured along
                axes specified by the Z syndrome graph

        Returns:
            logical_readout (int):
                logical readout value

            stabilizer_str (str):
                returns a string of the form
                "Z_{N}Z_{N-1}...Z_{0}"
        """
        readout_values = [int(q) for q in final_readout_string]
        readout_values = readout_values[::-1]  # [q_0,q_1,...,q_24]

        z_stabilizer = ""  # "Z_{N}Z_{N-1}..Z_{0}"

        for idx_list in self.geometry["mz"]:
            stabilizer_val = (
                sum(
                    [
                        readout_values[idx] if idx is not None else 0
                        for idx in idx_list[1:]
                    ]
                )  # [syn, left, right]
                % 2
            )
            z_stabilizer = str(stabilizer_val) + z_stabilizer

        logical_readout = readout_values[0]  # first row qubit

        return logical_readout, z_stabilizer

    def reset_x(self) -> None:
        """
        Initialize/reset to a logical |x+> state.
        Create a GHZ state: |+_L> := |0_L> + |1_L> = |00..0> + |11..1>
        """
        self.circ.reset(self.qregisters["data"])
        self.circ.h(self.qregisters["data"][0])
        for i in range(len(self.qregisters["data"]) - 1):
            self.circ.cx(self.qregisters["data"][i], self.qregisters["data"][i + 1])
        self.circ.barrier()

    def reset_z(self) -> None:
        """
        Initialize/reset to a logical |z+> state.
        """
        self.circ.reset(self.qregisters["data"])
        self.circ.barrier()

    def x(self) -> None:
        """
        Logical X operator on the topological qubit.
        Defined as the left-most column on the X Syndrome Graph.
        """
        self.circ.x(self.qregisters["data"])
        self.circ.barrier()

    def z(self) -> None:
        """
        Logical Z operator on the topological qubit.
        Defined as the top-most row on the Z Syndrome Graph.
        """
        self.circ.z(self.qregisters["data"][0])

    def x_c_if(self, classical: ClassicalRegister, val: int) -> None:
        """
        Classically conditioned logical X operator on the topological qubit.
        """
        self.circ.x(self.qregisters["data"]).c_if(classical, val)
        self.circ.barrier()

    def z_c_if(self, classical: ClassicalRegister, val: int) -> None:
        """
        Classically conditioned logical Z operator on the topological qubit.
        """
        self.circ.z(self.qregisters["data"][0]).c_if(classical, val)

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
            self.circ.cx(control, self.qregisters["data"])
            self.circ.barrier()
        elif target:
            self.circ.cx(self.qregisters["data"][0], target)

    def _readout_x_into_ancilla(self) -> None:
        """
        Convenience method to read-out the
        logical-X projection into an ancilla qubit.
        Uses the left-most column.
        """

        self.circ.reset(self.qregisters["ancilla"])

        # X Readout
        self.circ.h(self.qregisters["ancilla"])
        self.circ.cx(self.qregisters["ancilla"], self.qregisters["data"])
        self.circ.h(self.qregisters["ancilla"])

    def readout_x(self, readout_creg: Optional[ClassicalRegister] = None) -> None:
        """
        Convenience method to read-out the logical-X projection.
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

    def readout_z(self, readout_creg: Optional[ClassicalRegister] = None) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        Uses the top-most row (in this case just the first qubit).
        """
        if not readout_creg:
            self.params["num_readout"] += 1
            creg_name = self.name + "_readout_" + str(self.params["num_readout"])
            readout = ClassicalRegister(1, name=creg_name)

            self.circ.add_register(readout)

            self.cregisters[creg_name] = readout
            readout_creg = self.cregisters[creg_name]

        self.circ.reset(self.qregisters["ancilla"])
        self.circ.measure(self.qregisters["data"][0], readout_creg)
        self.circ.barrier()

    def lattice_readout_x(self) -> None:
        """
        Not applicable/relevant to the Rep Code.
        """
        self.params["num_lattice_readout"] += 1
        creg_name = (
            self.name + "_lattice_readout_" + str(self.params["num_lattice_readout"])
        )
        readout = ClassicalRegister(self.params["num_data"], name=creg_name,)
        self.circ.add_register(readout)
        self.cregisters[creg_name] = readout

        # measure along X
        self.circ.h(self.qregisters["data"])
        self.circ.measure(self.qregisters["data"], self.cregisters[creg_name])
        self.circ.barrier()

    def lattice_readout_z(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of Parity stabilizer measurments,
        as well as a logical Z readout.
        """
        self.params["num_lattice_readout"] += 1
        creg_name = (
            self.name + "_lattice_readout_" + str(self.params["num_lattice_readout"])
        )
        readout = ClassicalRegister(self.params["num_data"], name=creg_name,)

        self.circ.add_register(readout)

        self.cregisters[creg_name] = readout

        self.circ.measure(self.qregisters["data"], self.cregisters[creg_name])
        self.circ.barrier()

    def parse_readout(
        self, readout_string: str, readout_type: Optional[str] = None
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
        """
        Helper method to turn a result string (e.g. 0000 000 000) into an
        appropriate logical readout value and XOR-ed syndrome locations
        according to our grid coordinate convention.

        Args:
            readout_string (str):
                Readout like "0 000 000 000" (logical_readout syndrome_2 syndrome_1 syndrome_0)
                or of the form "0000 000 000" (lattice_readout syndrome_1 syndrome_0)
                A syndrome_i measurement "00..0" is of the form Z_{N}Z_{N-1}...Z_{0}

            readout_type (Optional[str]):
                "X" or "Z" needed to accurately parse a lattice readout to extract a final round of
                syndrome measurements and logical readout.
        Returns:
            logical_readout (int):
                logical readout value
            syndromes (Dict[str, List[TQubit]]]):
                key: syndrome type
                value: (time, row, col) of parsed syndrome hits (changes between consecutive rounds)
        """
        chunks = readout_string.split(" ")

        if len(chunks[0]) > 1:  # this is true when all data qubits are readout
            assert (
                readout_type == "Z"
            ), "Rep code currently only supports Z lattice readout."
            (
                logical_readout,
                final_stabilizer,
            ) = self.extract_final_stabilizer_and_logical_readout_z(chunks[0])
            chunks = [final_stabilizer,] + chunks[1:]
        else:
            logical_readout = int(chunks[0])
            chunks = chunks[1:]

        int_syndromes = [int(x, base=2) for x in chunks[::-1]]
        z_syndromes = [a ^ b for (a, b) in zip(int_syndromes, int_syndromes[1:])]

        Z = []
        for T, syndrome in enumerate(z_syndromes):
            for loc in range(int(self.params["num_syn"])):
                if syndrome & 1 << loc:
                    Z.append((float(T), 0.5 + loc, 0.0))
        return (
            logical_readout,
            {"Z": Z},
        )


class RepetitionQubit(TopologicalQubit):
    """
    A single logical repetition code qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend TopologicalQubit which extends QuantumCircuit.
    """

    lattice_type = _RepetitionLattice

    def stabilize(self) -> None:
        """
        Run a single round of stabilization (entangle and measure).
        """
        self.lattice.params["T"] += 1
        syndrome_readouts = ClassicalRegister(
            self.lattice.params["num_syn"],
            name=self.name + "_c{}".format(self.lattice.params["T"]),
        )
        self.lattice.cregisters[
            "syndrome{}".format(self.lattice.params["T"])
        ] = syndrome_readouts
        self.circ.add_register(syndrome_readouts)

        self.lattice.entangle()

        # measure syndromes
        self.circ.measure(
            self.lattice.qregisters["mz"], syndrome_readouts,
        )
        self.circ.reset(self.lattice.qregisters["mz"])
        self.circ.barrier()
