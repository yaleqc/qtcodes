"""
Repetition Code Encoder Classes
"""
from typing import Dict, List, Tuple, Optional, Type
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Qubit
from topological_codes.circuits.base import (
    _Stabilizer,
    LatticeError,
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


class _RepetitionLattice(_TopologicalLattice[TQubit]):
    """
    This class contains all the lattice geometry specifications regarding the Repetition Code.
    """

    def __init__(self, params: Dict[str, int], name: str, circ: QuantumCircuit):
        """
        Initializes this Topological Lattice class.

        Args:
            params (Dict[str,int]):
                Contains params such as d, where d is the number of
                physical "data" qubits lining a row or column of the lattice.
                Only odd d is possible!
            name (str):
                Useful when combining multiple TopologicalQubits together.
                Prepended to all registers.
            circ (QuantumCircuit):
                QuantumCircuit on top of which the topological qubit is built.
                This is often shared amongst multiple TQubits.
        """

        # validation
        required_params = ["d"]
        for required_param in required_params:
            if required_param not in params:
                raise LatticeError(f"Please include a {required_param} param.")

        # calculated params
        params["T"] = -1  # -1 until a stabilizer round is added!
        params["num_readout"] = -1  # -1 until a logical readout is performed!
        params["num_lattice_readout"] = -1  # -1 until a lattice readout is performed!
        params["num_data"] = params["d"]
        params["num_syn"] = params["d"] - 1

        # create registers
        qregisters: Dict[str, QuantumRegister] = {}  # quantum
        qregisters["data"] = QuantumRegister(params["num_data"], name=name + "_data")
        qregisters["mp"] = QuantumRegister(params["num_syn"], name=name + "_mp")

        cregisters: Dict[str, ClassicalRegister] = {}  # classical

        self.geometry: Dict[str, List[List[int]]] = {}

        super().__init__(circ, qregisters, cregisters, params, name)

    def _set_geometry(self):
        """
        Construct the lattice geometry for reuse across this class.

        Returns:
            geometry (Dict[str, List[List[int]]]):
                key: syndrome/plaquette type
                value: List of lists of qubit indices comprising one plaquette.
        """
        geometry = {"mp": []}

        for i in range(self.params["num_syn"]):
            syn = i
            left = i
            right = i + 1
            geometry["mp"].append([syn, left, right])

        self.geometry = geometry

    def gen_qubit_indices_and_stabilizers(
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
                "P_{N}P_{N-1}...P_{0}"
        """
        readout_values = [int(q) for q in final_readout_string]
        readout_values = readout_values[::-1]  # [q_0,q_1,...,q_24]

        p_stabilizer = ""  # "P_{N}P_{N-1}..P_{0}"

        for idx_list in self.geometry["mp"]:
            stabilizer_val = (
                sum(
                    [
                        readout_values[idx] if idx is not None else 0
                        for idx in idx_list[1:]
                    ]
                )  # [syn, left, right]
                % 2
            )
            p_stabilizer = str(stabilizer_val) + p_stabilizer

        logical_readout = 0
        for idx in range(self.params["d"]):
            logical_readout = (logical_readout + readout_values[idx]) % 2

        return logical_readout, p_stabilizer

    def parse_readout(
        self, readout_string: str, readout_type: str = "Z"
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
        """
        Helper method to turn a result string (e.g. 0000 000 000) into an
        appropriate logical readout value and XOR-ed syndrome locations
        according to our grid coordinate convention.

        Args:
            readout_string (str):
                Readout of the form "0000 000 000" (lattice_readout syndrome_1 syndrome_0)
        Returns:
            logical_readout (int):
                logical readout value
            syndromes (Dict[str, List[TQubit]]]):
                key: syndrome type
                value: (time, row, col) of parsed syndrome hits (changes between consecutive rounds)
        """
        assert readout_type == "Z", "Rep code only supports Z readout."

        chunks = readout_string.split(" ")

        (
            logical_readout,
            final_stabilizer,
        ) = self.extract_final_stabilizer_and_logical_readout_z(chunks[0])

        chunks = [final_stabilizer,] + chunks[1:]

        int_syndromes = [int(x, base=2) for x in chunks[::-1]]
        p_syndromes = [a ^ b for (a, b) in zip(int_syndromes, int_syndromes[1:])]

        P = []
        for T, syndrome in enumerate(p_syndromes):
            for loc in range(self.params["num_syn"]):
                if syndrome & 1 << loc:
                    P.append((float(T), 0.5 + loc // 2, 0.5 + loc % 2 * 2 - loc // 2))
        return (
            logical_readout,
            {"P": P},
        )


class RepetitionQubit(TopologicalQubit[TQubit]):
    """
    A single logical repetition code qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend TopologicalQubit which extends QuantumCircuit.
    """

    def __init__(
        self,
        params: Dict[str, int],
        name: str = "tqubit",
        circ: Optional[QuantumCircuit] = None,
    ) -> None:
        """
        Initializes this Topological Qubit class.

        Args:
            params (Dict[str,int]):
                Contains params such as d, where d is the number of
                physical "data" qubits within the rep code.
            name (str):
                Useful when combining multiple TopologicalQubits together.
                Prepended to all registers.
            circ (Optional[QuantumCircuit]):
                QuantumCircuit on top of which the topological qubit is built.
                This is often shared amongst multiple TQubits.
                If none is provided, then a new QuantumCircuit is initialized and stored.

        """
        # == None is necessary, as "not circ" is true for circ=QuantumCircuit()
        circ = QuantumCircuit() if circ is None else circ
        super().__init__(circ, name)
        self.lattice = _RepetitionLattice(params, name, circ)

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
            self.lattice.qregisters["mp"], syndrome_readouts,
        )
        self.circ.reset(self.lattice.qregisters["mp"])
        self.circ.barrier()

    def identity(self) -> None:
        """
        Inserts an identity on the data and syndrome qubits.
        This allows us to create an isolated noise model by inserting errors only on identity gates.
        """
        for register in self.lattice.qregisters.values():
            self.circ.id(register)
        self.circ.barrier()

    def identity_data(self) -> None:
        """
        Inserts an identity on the data qubits only.
        This allows us to create an isolated noise model by inserting errors only on identity gates.
        """
        self.circ.id(self.lattice.qregisters["data"])
        self.circ.barrier()

    def logical_x_plus_reset(self) -> None:
        """
        Initialize/reset to a logical |x+> state.
        """
        raise NotImplementedError("Not applicable/relevant to the Rep Code.")

    def logical_z_plus_reset(self) -> None:
        """
        Initialize/reset to a logical |z+> state.
        """
        self.circ.reset(self.lattice.qregisters["data"])
        self.circ.barrier()

    def logical_x(self) -> None:
        """
        Logical X operator on the topological qubit.
        Defined as the left-most column on the X Syndrome Graph.
        """
        self.circ.x(self.lattice.qregisters["data"])
        self.circ.barrier()

    def logical_z(self) -> None:
        """
        Logical Z operator on the topological qubit.
        Defined as the top-most row on the Z Syndrome Graph.
        """
        raise NotImplementedError("Not applicable/relevant to the Rep Code.")

    def readout_x(self) -> None:
        """
        Convenience method to read-out the logical-X projection.
        """
        raise NotImplementedError("Not applicable/relevant to the Rep Code.")

    def readout_z(self) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        """
        raise NotImplementedError("Not applicable/relevant to the Rep Code.")

    def lattice_readout_x(self) -> None:
        """
        Not applicable/relevant to the Rep Code.
        """
        raise NotImplementedError("Not applicable/relevant to the Rep Code.")

    def lattice_readout_z(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of Parity stabilizer measurments,
        as well as a logical Z readout.
        """
        self.lattice.params["num_lattice_readout"] += 1

        readout = ClassicalRegister(
            self.lattice.params["num_data"],
            name=self.name
            + "_lattice_readout_"
            + str(self.lattice.params["num_lattice_readout"]),
        )

        self.circ.add_register(readout)

        self.lattice.cregisters["lattice_readout"] = readout

        self.circ.measure(
            self.lattice.qregisters["data"], self.lattice.cregisters["lattice_readout"]
        )
        self.circ.barrier()

    def parse_readout(
        self, readout_string: str, readout_type: str = "Z"
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
        """
        Wrapper on helper method to turn a result string (e.g. 1 10100000 10010000) into an
        appropriate logical readout value and XOR-ed syndrome locations
        according to our grid coordinate convention.
        """
        return self.lattice.parse_readout(readout_string, readout_type)
