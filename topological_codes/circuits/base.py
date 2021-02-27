from abc import abstractmethod, ABCMeta
from typing import TypeVar, Tuple, Dict, List
from qiskit import QuantumRegister, QuantumCircuit

TQubit = TypeVar("TQubit")


class LatticeError(Exception):
    pass


class TopologicalLattice(metaclass=ABCMeta):
    @abstractmethod
    def __init__(
        self,
        d: int,
        data_register: QuantumRegister,
        mx_register: QuantumRegister,
        mz_register: QuantumRegister
    ) -> None:
        """
        Initializes an instance of the topological code lattice

        Args:
            d (int): surface code distance
            data_register (QuantumRegister): grouped register of all data qubits
            mx_register (QuantumRegister): grouped register of all measure-x qubits
            mz_register (QuantumRegister): grouped register of all measure-z qubits
        """
        pass

    @abstractmethod
    def entangle(self, circ: QuantumCircuit) -> None:
        pass

    @abstractmethod
    def parse_readout(
        self, readout_string: str
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
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

    @abstractmethod
    def __init__(
        self,
        data_register: QuantumRegister,
        mz_register: QuantumRegister,
        mx_register: QuantumRegister,
        ancilla: QuantumRegister
    ) -> None:
        super().__init__(data_register, mz_register, mx_register, ancilla)

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
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
        pass


# TODO: Implement TopologicalCircuit
class TopologicalCircuit:
    def __init__(self):
        pass