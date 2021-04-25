from topological_codes import (
    _Stabilizer,
    LatticeError,
    _TopologicalLattice,
    TopologicalQubit,
)
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, execute
from qiskit.circuit.quantumregister import Qubit
from typing import Dict, List, Tuple

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
    def __init__(self, circ: QuantumCircuit, params: Dict[str, int], name: str):
        # validation
        required_params = ["d"]
        for required_param in required_params:
            if required_param not in params:
                raise LatticeError(f"Please include a {required_param} param.")

        # calculated params
        params["T"] = -1  # -1 until a stabilizer round is added!
        params["num_data"] = params["d"]
        params["num_syn"] = params["d"] - 1

        # create registers
        qregisters: Dict[str, QuantumRegister] = {}  # quantum
        qregisters["data"] = QuantumRegister(params["num_data"], name=name + "_data")
        qregisters["mp"] = QuantumRegister(params["num_syn"], name=name + "_mp")

        cregisters: Dict[str, ClassicalRegister] = {}  # classical
        super().__init__(circ, qregisters, cregisters, params, name)

    def gen_qubit_indices_and_stabilizers(self):
        """
        Generates lattice blueprint for repetition code lattice with our
        chosen layout and numbering.
        """
        qubit_indices = []
        stabilizers = []

        for i in range(self.params["num_syn"]):
            syn = self.qregisters["mp"][i]
            left = self.qregisters["data"][i]
            right = self.qregisters["data"][i + 1]
            qubit_indices.append([syn, left, right])
            stabilizers.append(_Parity)
        return qubit_indices, stabilizers

    def parse_readout(self, readout_string: str):
        pass


class RepetitionQubit(TopologicalQubit):
    """
    A single logical repetition code qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend TopologicalQubit which extends QuantumCircuit.
    """

    def __init__(
        self, params: Dict[str, int], name: str = "tqubit", circ: QuantumCircuit = None,
    ) -> None:
        """
        Initializes a new QuantumCircuit for this logical qubit and calculates
        the underlying surface code lattice ordering.
        
        Args:
            d (int): Number of physical "data" qubits. Only odd d is possible!
        """
        circ = circ if circ else QuantumCircuit()
        super().__init__(circ, name)
        self.lattice = _RepetitionLattice(circ, params, name)

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
        [
            self.circ.id(x)
            for register in self.lattice.qregisters.values()
            for x in register
        ]
        self.circ.barrier()

    def identity_data(self) -> None:
        [self.circ.id(x) for x in self.lattice.qregisters["data"]]
        self.circ.barrier()

    def hadamard_reset(self) -> None:
        raise NotImplementedError("This has not been implemented yet.")

    def logical_x(self) -> None:
        self.circ.x(self.lattice.qregisters["data"])
        self.circ.barrier()

    def logical_z(self) -> None:
        raise NotImplementedError("This has not been implemented yet.")

    def readout_z(self) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        """
        readout = ClassicalRegister(
            self.lattice.params["num_data"], name=self.name + "_readout"
        )

        # try adding readout cregister
        # this will throw an error if a "readout" register is already a part of the circ
        # TODO: add functionality to have multiple readout registers
        self.circ.add_register(readout)
        self.lattice.cregisters["readout"] = readout
        self.circ.measure(
            self.lattice.qregisters["data"], self.lattice.cregisters["readout"]
        )
        self.circ.barrier()

    def readout_x(self) -> None:
        raise NotImplementedError("This has not been implemented yet.")

    def parse_readout(
        self, readout_string: str
    ):  # TODO: -> Tuple[int, Dict[str, List[TQubit]]]:
        return self.lattice.parse_readout(readout_string)

