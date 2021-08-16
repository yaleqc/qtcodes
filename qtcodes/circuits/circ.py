"""
Topological Circuit and Register
"""
from typing import Union, Dict, cast, Optional, Any, Tuple, List, Type
from qiskit import QuantumCircuit, QuantumRegister
from qiskit.circuit.classicalregister import ClassicalRegister
from qtcodes.circuits.base import TopologicalQubit
from qtcodes.circuits.constants import *


class TopologicalRegister:
    """
    A blueprint for a TopologicalRegister that stores topological qubit(s)
    """

    def __init__(
        self,
        num_tqubits: int,
        circ: QuantumCircuit = None,
        ctype: str = XXZZ,
        params: Optional[Dict[str, int]] = None,
        name: str = "treg",
    ):
        """
        Args:
            num_tqubits (int):
                The number of topological qubits to be stored in the register.
            circ (QuantumCircuit):
                QuantumCircuit on top of which the topological qubit is built.
                This is often shared amongst multiple TQubits.
                If none is provided, then a new QuantumCircuit is initialized and stored.
            type (str):
                Specifies the type of TopologicalQubit
                Must be a value of the blueprint dictionary
            params (Dict[str,int]):
                 Contains params such as d, where d is the number of
                 physical "data" qubits lining a row or column of the lattice.
            name (str):
                Useful when combining multiple TopologicalQubits together.
                Prepended to all registers.

        """
        self.params = params if params else {"d": 3}
        self.name = name

        self.n = 0

        # == None is necessary, as "not circ" is true for circ=QuantumCircuit()
        self.circ = QuantumCircuit() if circ is None else circ

        if ctype not in blueprint:
            raise ValueError(
                "Please choose a Topological Qubit type from: "
                + str(list(blueprint.keys()))
            )
        self.tqubit_type = blueprint[ctype]
        self.tqubits: Dict[str, Dict[int, Type[Any]]] = {}
        for _ in range(num_tqubits):
            self.add_tqubit("data")

    def add_tqubit(self, sub_register: str) -> None:
        if sub_register not in self.tqubits:
            self.tqubits[sub_register] = {}
        self.tqubits[sub_register][self.n] = self.tqubit_type(
            params=self.params, name=self.name + "_" + str(self.n), circ=self.circ
        )
        self.n += 1

    def add_tqubits(self, sub_register: str, num_tqubits: int) -> None:
        for _ in range(num_tqubits):
            self.add_tqubit(sub_register)

    def __getitem__(self, key: Union[str, int]):
        """
        Allows us to return the nth element of TopologicalRegister as a list.
        """
        if isinstance(key, str):
            key = str(key)
            return self.tqubits[key]
        elif isinstance(key, int):
            key = int(key)
            for _, sub_reg in self.tqubits.items():
                if key in sub_reg:
                    return sub_reg[key]
        raise ValueError("Key not found!")


class TopologicalCircuit:
    """
    TopologicalCircuit is like a QuantumCircuit built on Topological Qubits.
    Shares the same QuantumCircuit object created in TopologicalRegister.
    """

    def __init__(self, treg: TopologicalRegister):
        self.treg = treg
        self.qreg: Dict[str, QuantumRegister] = {}
        self.creg: Dict[str, ClassicalRegister] = {}
        self.circ = treg.circ

    def add_creg(self, size=None, name=None, bits=None, override: bool = False) -> None:
        if name in self.creg and not override:
            return
        creg = ClassicalRegister(size=size, name=name, bits=bits)
        self.creg[name] = creg
        self.circ.add_register(creg)

    def add_qreg(self, size=None, name=None, bits=None, override: bool = False) -> None:
        if name in self.qreg and not override:
            return
        qreg = QuantumRegister(size=size, name=name, bits=bits)
        self.qreg[name] = qreg
        self.circ.add_register(qreg)

    def _get_index(self, tqubit: Union[TopologicalQubit, int]) -> TopologicalQubit:
        """
        Takes in either a TopologicalQubit or an int, and returns a TopologicalQubit.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg

        Returns:
            tqubit (TopologicalQubit):
                Returns the corresponding TopologicalQubit from treg
        """
        if isinstance(tqubit, int):
            tqubit = cast(int, tqubit)
            tqubit = self.treg[tqubit]
        tqubit = cast(TopologicalQubit, tqubit)
        return tqubit

    def stabilize(self, tqubit: Union[TopologicalQubit, int]):
        """
        Run a single round of stabilization (entangle and measure) on the tqubit.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.stabilize()

    def id(self, tqubit: Union[TopologicalQubit, int]) -> None:
        """
        Inserts an identity on the data and syndrome qubits.
        This allows us to create an isolated noise model by inserting errors only on identity gates.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.id()

    def id_data(self, tqubit: Union[TopologicalQubit, int]) -> None:
        """
        Inserts an identity on the data qubits only.
        This allows us to create an isolated noise model by inserting errors only on identity gates.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.id_data()

    def reset_x(self, tqubit: Union[TopologicalQubit, int]) -> None:
        """
        Initialize/reset to a logical |x+> state on the tqubit.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.reset_x()

    def reset_z(self, tqubit: Union[TopologicalQubit, int]):
        """
        Initialize/reset to a logical |z+> state on the tqubit.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.reset_z()

    def x(self, tqubit: Union[TopologicalQubit, int]):
        """
        Logical X operator on the topological qubit.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.x()

    def z(self, tqubit: Union[TopologicalQubit, int]):
        """
        Logical Z operator on the topological qubit.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.z()

    def cx(
        self,
        control: Union[TopologicalQubit, int],
        target: Union[TopologicalQubit, int],
    ):

        # get qubits
        control = self._get_index(control)
        target = self._get_index(target)
        if "ancilla" not in self.treg.tqubits:
            self.treg.add_tqubit("ancilla")
        ancilla = cast(TopologicalQubit, list(self.treg["ancilla"].values())[-1])

        # prepare bits
        self.add_creg(1, "m1")
        self.add_creg(1, "m2")
        self.add_creg(1, "m3")

        # prepare ancilla
        ancilla.reset_x()

        # Z (control) x Z (ancilla)
        self.add_qreg(1, "cnot_readout")

        readout = self.qreg["cnot_readout"][0]

        # Z x Z
        self.circ.reset(readout)
        control.cx(target=readout)
        ancilla.cx(target=readout)
        self.circ.measure(readout, self.creg["m1"][0])

        # X x X
        self.circ.reset(readout)
        self.circ.h(readout)
        target.cx(control=readout)
        ancilla.cx(control=readout)
        self.circ.h(readout)
        self.circ.measure(readout, self.creg["m2"][0])

        # Z
        ancilla.readout_z(readout_creg=self.creg["m3"])

        # classical conditioned
        control.z_c_if(self.creg["m2"], 1)
        target.x_c_if(self.creg["m1"], 1)
        target.x_c_if(self.creg["m3"], 1)

    def measure_x(
        self,
        tqubit: Union[TopologicalQubit, int],
        readout_creg: Optional[ClassicalRegister] = None,
    ):
        """
        Convenience method to read-out the logical-X projection.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.readout_x(readout_creg=readout_creg)

    def measure_z(
        self,
        tqubit: Union[TopologicalQubit, int],
        readout_creg: Optional[ClassicalRegister] = None,
    ):
        """
        Convenience method to read-out the logical-Z projection.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.readout_z(readout_creg=readout_creg)

    def measure_lattice_x(self, tqubit: Union[TopologicalQubit, int]):
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of stabilizer measurments,
        as well as a logical X readout.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.lattice_readout_x()

    def measure_lattice_z(self, tqubit: Union[TopologicalQubit, int]):
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of stabilizer measurments,
        as well as a logical Z readout.

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
        """
        tqubit = self._get_index(tqubit)
        tqubit.lattice_readout_z()

    def parse_readout(
        self,
        tqubit: Union[TopologicalQubit, int],
        readout_string: str,
        readout_type: Optional[str] = "Z",
    ) -> Tuple[int, Dict[str, List[Any]]]:
        """
        Helper method to turn a result string (e.g. 1 10100000 10010000) into an
        appropriate logical readout value and XOR-ed syndrome locations
        according to our grid coordinate convention, based on the topological qubit of choice.

        The implementation varies with different topological qubits,
        but here's an example from the rotated surface code:

        Args:
            tqubit (Union[TopologicalQubit, int]):
                Either already a TopologicalQubit or an int index in treg
            readout_string (str):
                Readout of the form "0 00000000 00000000" (logical_readout syndrome_1 syndrome_0)
                or of the form "000000000 00000000 00000000" (lattice_readout syndrome_1 syndrome_0)
        Returns:
            logical_readout (int):
                logical readout value
            syndromes (Dict[str, List[TQubit]]]):
                key: syndrome type
                value: (time, row, col) of parsed syndrome hits (changes between consecutive rounds)
        """
        tqubit = self._get_index(tqubit)
        return tqubit.parse_readout(readout_string, readout_type)

    def draw(self, **kwargs):
        """
        Convenience method to draw underlying quantum circuit.
        """
        return self.circ.draw(**kwargs)

    def __str__(self):
        return self.circ.__str__()
