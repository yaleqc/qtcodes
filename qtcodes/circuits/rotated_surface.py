"""
Rotated Surface Code Encoder Classes
"""
from abc import abstractmethod, ABCMeta
from typing import Dict, List, Tuple, Optional, Type
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister
from qiskit.circuit import Qubit

from qtcodes.circuits.base import (
    LatticeError,
    _TopologicalLattice,
    TopologicalQubit,
    _Stabilizer,
)

TQubit = Tuple[float, float, float]


class _RotatedLattice(_TopologicalLattice[TQubit], metaclass=ABCMeta):
    """
    This class contains all the lattice geometry specifications
    regarding the XXZZ (CSS) Rotated Surface Code.
    """

    @property
    @abstractmethod
    def stabilizer_shortnames(self) -> Dict[str, Type[_Stabilizer]]:
        """
        key: "mx", "mz"
        values: _Stabilizer subclass
        """

    def __init__(self, params: Dict[str, float], name: str, circ: QuantumCircuit):
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
        self.geometry: Dict[str, List[List[Optional[int]]]] = {}
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

        # validation
        if self.params["d"] % 2 != 1:
            raise LatticeError("Surface code distance must be odd!")

        # calculated params
        self.params["T"] = -1  # -1 until a stabilizer round is added!
        self.params["num_readout"] = -1  # -1 until a logical readout is performed!
        self.params[
            "num_lattice_readout"
        ] = -1  # -1 until a lattice readout is performed!
        self.params["num_data"] = self.params["d"] ** 2
        self.params["num_syn"] = (self.params["d"] ** 2 - 1) // 2

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
            self.params["num_syn"], name=self.name + "_mz"
        )
        self.qregisters["mx"] = QuantumRegister(
            self.params["num_syn"], name=self.name + "_mx"
        )
        self.qregisters["ancilla"] = QuantumRegister(1, name=self.name + "_ancilla")

    def _set_geometry(self) -> None:
        """
        Construct the lattice geometry for reuse across this class.

        Returns:
            geometry (Dict[str, List[List[int]]]):
                key: syndrome/plaquette type
                value: List of lists of qubit indices comprising one plaquette.
        """
        geometry: Dict[str, List[List[Optional[int]]]] = {"mx": [], "mz": []}
        d = int(self.params["d"])
        per_row_x = (d - 1) // 2
        per_row_z = (d + 1) // 2
        # mx geometry

        # good typing
        top_l: Optional[int] = None
        top_r: Optional[int] = None
        bot_l: Optional[int] = None
        bot_r: Optional[int] = None

        for syn in range(int(self.params["num_syn"])):
            row = syn // per_row_x
            offset = syn % per_row_x
            start = (row - 1) * d
            row_parity = row % 2

            if row == 0:  # First row
                top_l, top_r = None, None
                bot_l = syn * 2
                bot_r = syn * 2 + 1
            elif row == d:  # Last row
                bot_l, bot_r = None, None
                top_l = syn * 2 + 1
                top_r = syn * 2 + 2
            else:
                top_l = start + (offset * 2) + row_parity
                top_r = start + (offset * 2) + row_parity + 1
                bot_l = start + d + (offset * 2) + row_parity
                bot_r = start + d + (offset * 2) + row_parity + 1

            geometry["mx"].append([syn, top_l, top_r, bot_l, bot_r])

        for syn in range(int(self.params["num_syn"])):
            row = syn // per_row_z
            offset = syn % per_row_z
            start = row * d
            row_parity = row % 2

            top_l = start + (offset * 2) - row_parity
            top_r = start + (offset * 2) - row_parity + 1
            bot_l = start + d + (offset * 2) - row_parity
            bot_r = start + d + (offset * 2) - row_parity + 1

            # Overwrite edge column syndromes
            if row_parity == 0 and offset == per_row_z - 1:  # Last column
                top_r, bot_r = None, None
            elif row_parity == 1 and offset == 0:  # First column
                top_l, bot_l = None, None

            geometry["mz"].append([syn, top_l, top_r, bot_l, bot_r])
        self.geometry = geometry

    def _gen_qubit_indices_and_stabilizers(
        self,
    ) -> Tuple[List[List[Qubit]], List[Type[_Stabilizer]]]:
        """
        Generates lattice blueprint for rotated surface code lattice with our
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
            stabilizer_cls = self.stabilizer_shortnames[stabilizer]
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

    def entangle_x(self):
        """
        Build/entangle just the X Syndrome circuit.
        """
        num_x = len(self.qregisters["mx"])
        self.entangle(self.qubit_indices[:num_x], self.stabilizers[:num_x])

    def entangle_z(self):
        """
        Build/entangle just the Z Syndrome circuit.
        """
        num_x = len(self.qregisters["mx"])
        self.entangle(self.qubit_indices[num_x:], self.stabilizers[num_x:])

    def extract_final_stabilizer_and_logical_readout_x(
        self, final_readout_string: str, previous_syndrome_string: str
    ) -> Tuple[int, str]:
        """
        Extract final X syndrome measurements and logical X readout from
        lattice readout along the X syndrome graph.

        Args:
            final_readout_string (str):
                readout string of length equal to the number of data qubits
                contains the readout values of each data qubit measured along
                axes specified by the X syndrome graph

            previous_syndrome_string (str):
                syndrome readout string form the previous round
                of syndrome/stabilizer measurements


        Returns:
            logical_readout (int):
                logical readout value

            stabilizer_str (str):
                returns a string of the form
                "X_{N}X_{N-1}...X_{0}Z_{N}Z_{N-1}...Z_{0}",
                where Z_{N}Z_{N-1}...Z_{0} is copied from the previous Z syndrome
                readout stored in previous_syndrome_string
        """
        readout_values = [int(q) for q in final_readout_string]
        readout_values = readout_values[::-1]  # [q_0,q_1,...,q_24]

        x_stabilizer = ""  # "X_{N}X_{N-1}..X_{0}"

        for idx_list in self.geometry["mx"]:
            stabilizer_val = (
                sum(
                    [
                        readout_values[idx] if idx is not None else 0
                        for idx in idx_list[1:]
                    ]
                )  # [mx, top_l, top_r, bot_l, bot_r]
                % 2
            )
            x_stabilizer = str(stabilizer_val) + x_stabilizer

        stabilizer_str = (
            x_stabilizer + previous_syndrome_string[int(self.params["num_syn"]) :]
        )
        # X_{N}X_{N-1}...X_{0}Z_{N}Z_{N-1}...Z_{0}, where
        # Z_{N}Z_{N-1}...Z_{0} is copied from previous syndrome measurement string

        logical_readout = 0
        for idx in range(0, int(self.params["num_data"]), int(self.params["d"])):
            logical_readout = (logical_readout + readout_values[idx]) % 2

        return logical_readout, stabilizer_str

    def extract_final_stabilizer_and_logical_readout_z(
        self, final_readout_string: str, previous_syndrome_string: str
    ) -> Tuple[int, str]:
        """
        Extract final Z syndrome measurements and logical Z readout from
        lattice readout along the Z syndrome graph.

        Args:
            final_readout_string (str):
                readout string of length equal to the number of data qubits
                contains the readout values of each data qubit measured along
                axes specified by the Z syndrome graph

            previous_syndrome_string (str):
                syndrome readout string form the previous round
                of syndrome/stabilizer measurements


        Returns:
            logical_readout (int):
                logical readout value

            stabilizer_str (str):
                returns a string of the form
                "X_{N}X_{N-1}...X_{0}Z_{N}Z_{N-1}...Z_{0}",
                where X_{N}X_{N-1}...X_{0} is copied from the previous X syndrome
                readout stored in previous_syndrome_string
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
                )  # [mx, top_l, top_r, bot_l, bot_r]
                % 2
            )
            z_stabilizer = str(stabilizer_val) + z_stabilizer

        stabilizer_str = (
            previous_syndrome_string[: int(self.params["num_syn"])] + z_stabilizer
        )
        # X_{N}X_{N-1}...X_{0}Z_{N}Z_{N-1}...Z_{0}, where
        # X_{N}X_{N-1}...X_{0} is copied from previous syndrome measurement string

        logical_readout = 0
        for idx in range(int(self.params["d"])):
            logical_readout = (logical_readout + readout_values[idx]) % 2

        return logical_readout, stabilizer_str

    def extract_final_stabilizer_and_logical_readout(
        self,
        final_readout_string: str,
        previous_syndrome_string: str,
        readout_type: str,
    ):
        """
        Wrapper on extract_final_stabilizer_and_logical_readout_x/z.
        """
        if readout_type == "X":
            return self.extract_final_stabilizer_and_logical_readout_x(
                final_readout_string, previous_syndrome_string
            )
        elif readout_type == "Z":
            return self.extract_final_stabilizer_and_logical_readout_z(
                final_readout_string, previous_syndrome_string
            )

    @abstractmethod
    def reset_x(self) -> None:
        """
        Initialize/reset to a logical |x+> state.
        """

    @abstractmethod
    def reset_z(self) -> None:
        """
        Initialize/reset to a logical |z+> state.
        """

    @abstractmethod
    def x(self) -> None:
        """
        Logical X operator on the qubit.
        Uses the left-most column.
        """

    @abstractmethod
    def z(self) -> None:
        """
        Logical Z operator on the qubit.
        Uses the top-most row.
        """

    @abstractmethod
    def readout_x(self, readout_creg: Optional[ClassicalRegister] = None) -> None:
        """
        Convenience method to read-out the logical-X projection.
        Uses the left-most column.
        """

    @abstractmethod
    def readout_z(self, readout_creg: Optional[ClassicalRegister] = None) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        Uses the top-most row.
        """

    @abstractmethod
    def lattice_readout_x(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of stabilizer measurments,
        as well as a logical X readout.
        """

    @abstractmethod
    def lattice_readout_z(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of stabilizer measurments,
        as well as a logical Z readout.
        """

    def parse_readout(
        self, readout_string: str, readout_type: Optional[str] = None
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
        """
        Helper method to turn a result string (e.g. 1 10100000 10010000) into an
        appropriate logical readout value and XOR-ed syndrome locations
        according to our grid coordinate convention.

        Args:
            readout_string (str):
                Readout of the form "0 00000000 00000000" (logical_readout syndrome_1 syndrome_0)
                or of the form "000000000 00000000 00000000" (lattice_readout syndrome_1 syndrome_0)

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
        d = int(self.params["d"])
        num_syn = int(self.params["num_syn"])

        if len(chunks[0]) > 1:  # this is true when all data qubits are readout
            assert readout_type is not None
            (
                logical_readout,
                final_stabilizer,
            ) = self.extract_final_stabilizer_and_logical_readout(
                chunks[0], chunks[1], readout_type
            )
            chunks = [final_stabilizer,] + chunks[1:]
        else:
            logical_readout = int(chunks[0])
            chunks = chunks[1:]

        int_syndromes = [int(x, base=2) for x in chunks[::-1]]
        xor_syndromes = [a ^ b for (a, b) in zip(int_syndromes, int_syndromes[1:])]

        mask_z = "1" * num_syn
        mask_x = mask_z + "0" * num_syn
        x_syndromes = [(x & int(mask_x, base=2)) >> num_syn for x in xor_syndromes]
        z_syndromes = [x & int(mask_z, base=2) for x in xor_syndromes]

        X = []
        per_row_x = d // 2
        for T, syndrome in enumerate(x_syndromes):
            for loc in range(num_syn):
                if syndrome & 1 << loc:
                    row = -0.5 + loc // per_row_x
                    col = (0.5 + (loc // per_row_x) % 2) + (loc % per_row_x) * 2
                    X.append((float(T), row, col))

        Z = []
        per_row_z = d // 2 + 1
        for T, syndrome in enumerate(z_syndromes):
            for loc in range(num_syn):
                if syndrome & 1 << loc:
                    row = 0.5 + loc // per_row_z
                    col = (0.5 - (loc // per_row_z) % 2) + (loc % per_row_z) * 2
                    Z.append((float(T), row, col))

        return (
            logical_readout,
            {"X": X, "Z": Z,},
        )


class RotatedQubit(TopologicalQubit[TQubit], metaclass=ABCMeta):
    """
    A single logical surface code qubit.
    """

    def stabilize(self) -> None:
        """
        Run a single round of stabilization (entangle and measure).
        """
        self.lattice.params["T"] += 1
        syndrome_readouts = ClassicalRegister(
            self.lattice.params["num_syn"] * 2,
            name=self.name + "_c{}".format(self.lattice.params["T"]),
        )
        self.lattice.cregisters[
            "syndrome{}".format(self.lattice.params["T"])
        ] = syndrome_readouts
        self.circ.add_register(syndrome_readouts)

        self.lattice.entangle()

        # measure syndromes
        self.circ.measure(
            self.lattice.qregisters["mz"],
            syndrome_readouts[0 : self.lattice.params["num_syn"]],
        )
        self.circ.measure(
            self.lattice.qregisters["mx"],
            syndrome_readouts[
                self.lattice.params["num_syn"] : self.lattice.params["num_syn"] * 2
            ],
        )
        self.circ.reset(self.lattice.qregisters["mz"])
        self.circ.reset(self.lattice.qregisters["mx"])
        self.circ.barrier()
