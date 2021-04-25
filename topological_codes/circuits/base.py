from abc import abstractmethod, ABCMeta
from typing import TypeVar, Tuple, Dict, List, Generic
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit.quantumregister import Qubit

TQubit = TypeVar("TQubit")


class LatticeError(Exception):
    pass


class _Stabilizer(metaclass=ABCMeta):
    """
    A blueprint for stabilizer classes, such as plaquettes for surface codes.
    """

    def __init__(self, circ: QuantumCircuit, qubit_indices: List[List[Qubit]]):
        self.circ = circ
        self.qubit_indices = qubit_indices

    @abstractmethod
    def entangle(self):
        pass


class _TopologicalLattice(metaclass=ABCMeta):
    def __init__(
        self,
        circ: QuantumCircuit,
        qregisters: Dict[str, QuantumRegister],
        cregisters: Dict[str, ClassicalRegister],
        params: Dict[str, int],
        name: str,
    ):
        self.circ = circ
        self.qregisters = qregisters
        self.cregisters = cregisters
        self.name = name

        # add registerse to circ
        registers = list(self.qregisters.values()) + list(self.cregisters.values())
        self.circ.add_register(*registers)

        self.params = params

        self.qubit_indices, self.stabilizers = self.gen_qubit_indices_and_stabilizers()

    @abstractmethod
    def gen_qubit_indices_and_stabilizers(
        self,
    ):  # TODO: adding -> Tuple[List[List[Qubit]], List[_Stabilizer]] is giving a typing error in self.entangle
        pass

    def entangle(self, qubit_indices=None, stabilizers=None) -> None:
        qubit_indices = qubit_indices if qubit_indices else self.qubit_indices
        stabilizers = stabilizers if stabilizers else self.stabilizers

        for i in range(len(stabilizers)):
            stabilizer = stabilizers[i](self.circ, qubit_indices[i])
            stabilizer.entangle()
            self.circ.barrier()

    @abstractmethod
    def parse_readout(
        self, readout_string: str
    ):  # TODO: -> Tuple[int, Dict[str, List[TQubit]]]:
        """
        Helper method to turn a result string (e.g. 1 10100000 10010000) into an
        appropriate logical readout value and XOR-ed syndrome locations
        according to our grid coordinate convention.
        """
        pass


class TopologicalQubit(QuantumCircuit, metaclass=ABCMeta):
    """
    A single topological code logical qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend QuantumCircuit.
    """

    def __init__(self, name: str):
        super().__init__()
        self.name = name

    @abstractmethod
    def stabilize(self) -> None:
        pass

    @abstractmethod
    def identity(self) -> None:
        pass

    @abstractmethod
    def identity_data(self) -> None:
        pass

    @abstractmethod
    def hadamard_reset(self) -> None:
        pass

    @abstractmethod
    def logical_x(self) -> None:
        pass

    @abstractmethod
    def logical_z(self) -> None:
        pass

    @abstractmethod
    def readout_z(self) -> None:
        pass

    @abstractmethod
    def readout_x(self) -> None:
        pass

    @abstractmethod
    def parse_readout(
        self, readout_string: str
    ):  # TODO: -> Tuple[int, Dict[str, List[TQubit]]]:
        pass
