from topological_codes import (
    _Stabilizer,
    LatticeError,
    _TopologicalLattice,
    TopologicalQubit,
)
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, execute
from qiskit.circuit.quantumregister import Qubit
from typing import Dict, List, Tuple, Optional

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


class _XZZXLattice(_TopologicalLattice[TQubit]):
    """
    This class contains all the lattice geometry specifications regarding the XZZX Rotated Surface Code.
    """

    def __init__(self, params: Dict[str, int], name: str, circ: QuantumCircuit):
        """
        Initializes this Topological Lattice class.
        
        Args:
            params (Dict[str,int]): 
                Contains params such as d, where d is the number of physical "data" qubits lining a row or column of the lattice. 
                Only odd d is possible!
            name (str):
                Useful when combining multiple TopologicalQubits together. Prepended to all registers.
            circ (QuantumCircuit):
                QuantumCircuit on top of which the topological qubit is built. This is often shared amongst multiple TQubits.
        """

        # validation
        required_params = ["d"]
        for required_param in required_params:
            if required_param not in params:
                raise LatticeError(f"Please include a {required_param} param.")
        if params["d"] % 2 != 1:
            raise LatticeError("Surface code distance must be odd!")

        # calculated params
        params["T"] = -1  # -1 until a stabilizer round is added!
        params["num_readout"] = -1  # -1 until a logical readout is performed!
        params["num_lattice_readout"] = -1  # -1 until a lattice readout is performed!
        params["num_data"] = params["d"] ** 2
        params["num_syn"] = (params["d"] ** 2 - 1) // 2

        # create registers
        qregisters: Dict[str, QuantumRegister] = {}  # quantum
        qregisters["data"] = QuantumRegister(params["num_data"], name=name + "_data")
        qregisters["mb"] = QuantumRegister(params["num_syn"], name=name + "_mb")
        qregisters["ma"] = QuantumRegister(params["num_syn"], name=name + "_ma")
        qregisters["ancilla"] = QuantumRegister(1, name=name + "_ancilla")

        cregisters: Dict[str, ClassicalRegister] = {}  # classical
        super().__init__(circ, qregisters, cregisters, params, name)

    def _geometry(self):
        """
        Construct the lattice geometry for reuse across this class.

        Returns:
            geometry (Dict[str, List[List[int]]]): 
                key: syndrome/plaquette type
                value: List of lists of qubit indices comprising one plaquette.
        """
        geometry = {"ma": [], "mb": []}

        d = self.params["d"]
        per_row_x = (d - 1) // 2
        per_row_z = (d + 1) // 2
        # mx geometry
        for syn in range(self.params["num_syn"]):
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

            geometry["ma"].append([syn, top_l, top_r, bot_l, bot_r])

        for syn in range(self.params["num_syn"]):
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

            geometry["mb"].append([syn, top_l, top_r, bot_l, bot_r])
        return geometry

    def gen_qubit_indices_and_stabilizers(self):
        """
        Generates lattice blueprint for rotated surface code lattice with our
        chosen layout and numbering.

        Returns:
            qubit_indices (List[List[Qubit]]):
                List of lists of Qubits that comprise each plaquette.

            stabilizers (List[_Stabilizer]):
                List of stabilizers for each plaquette.
        """
        self.geometry = self._geometry()

        qubit_indices = []
        stabilizers = []
        for stabilizer, idx_lists in self.geometry.items():
            stabilizer_cls = _XZZX
            for idx_list in idx_lists:
                syn = self.qregisters[stabilizer][idx_list[0]]
                plaquette = [
                    self.qregisters["data"][idx] if idx != None else None
                    for idx in idx_list[1:]
                ]
                plaquette = [syn,] + plaquette
                qubit_indices.append(plaquette)
                stabilizers.append(stabilizer_cls)
        return qubit_indices, stabilizers

    def entangle_a(self):
        """
        Build/entangle just the A Syndrome circuit.
        """
        num_x = len(self.qregisters["ma"])
        self.entangle(self.qubit_indices[:num_x], self.stabilizers[:num_x])

    def entangle_b(self):
        """
        Build/entangle just the B Syndrome circuit.
        """
        num_x = len(self.qregisters["ma"])
        self.entangle(self.qubit_indices[num_x:], self.stabilizers[num_x:])

    def extract_final_stabilizer_and_logical_readout_x(
        self, final_readout_string: str, previous_syndrome_string: str
    ) -> Tuple[int, str]:
        """
        Extract final A syndrome measurements and logical X readout from 
        lattice readout along the A syndrome graph.

        Args:
            final_readout_string (str):
                readout string of length equal to the number of data qubits
                contains the readout values of each data qubit measured along
                axes specified by the A syndrome graph
            
            previous_syndrome_string (str):
                syndrome readout string form the previous round 
                of syndrome/stabilizer measurements


        Returns:
            logical_readout (int):
                logical readout value
            
            stabilizer_str (str):
                returns a string of the form 
                "A_{N}A_{N-1}...A_{0}B_{N}B_{N-1}...B_{0}", 
                where B_{N}B_{N-1}...B_{0} is copied from the previous B syndrome 
                readout stored in previous_syndrome_string
        """
        readout_values = [int(q) for q in final_readout_string]
        readout_values = readout_values[::-1]  # [q_0,q_1,...,q_24]

        a_stabilizer = ""  # "A_{N}A_{N-1}..A_{0}"

        for idx_list in self.geometry["ma"]:
            stabilizer_val = (
                sum(
                    [readout_values[idx] if idx != None else 0 for idx in idx_list[1:]]
                )  # [mx, top_l, top_r, bot_l, bot_r]
                % 2
            )
            a_stabilizer = str(stabilizer_val) + a_stabilizer

        stabilizer_str = (
            a_stabilizer + previous_syndrome_string[self.params["num_syn"] :]
        )  # A_{N}A_{N-1}...A_{0}B_{N}B_{N-1}...B_{0}, where B_{N}B_{N-1}...B_{0} is copied from previous syndrome measurement string

        logical_readout = 0
        for idx in range(0, self.params["num_data"], self.params["d"]):
            logical_readout = (logical_readout + readout_values[idx]) % 2

        return logical_readout, stabilizer_str

    def extract_final_stabilizer_and_logical_readout_z(
        self, final_readout_string: str, previous_syndrome_string: str
    ) -> Tuple[int, str]:
        """
        Extract final B syndrome measurements and logical Z readout from 
        lattice readout along the B syndrome graph.

        Args:
            final_readout_string (str):
                readout string of length equal to the number of data qubits
                contains the readout values of each data qubit measured along
                axes specified by the B syndrome graph
            
            previous_syndrome_string (str):
                syndrome readout string form the previous round 
                of syndrome/stabilizer measurements


        Returns:
            logical_readout (int):
                logical readout value
            
            stabilizer_str (str):
                returns a string of the form 
                "A_{N}A_{N-1}...A_{0}B_{N}B_{N-1}...B_{0}", 
                where A_{N}A_{N-1}...A_{0} is copied from the previous A syndrome 
                readout stored in previous_syndrome_string
        """
        readout_values = [int(q) for q in final_readout_string]
        readout_values = readout_values[::-1]  # [q_0,q_1,...,q_24]

        b_stabilizer = ""  # "B_{N}B_{N-1}..B_{0}"

        for idx_list in self.geometry["mb"]:
            stabilizer_val = (
                sum(
                    [readout_values[idx] if idx != None else 0 for idx in idx_list[1:]]
                )  # [mx, top_l, top_r, bot_l, bot_r]
                % 2
            )
            b_stabilizer = str(stabilizer_val) + b_stabilizer

        stabilizer_str = (
            previous_syndrome_string[: self.params["num_syn"]] + b_stabilizer
        )  # A_{N}A_{N-1}...A_{0}B_{N}B_{N-1}...B_{0}, where A_{N}A_{N-1}...A_{0} is copied from previous syndrome measurement string

        logical_readout = 0
        for id in range(self.params["d"]):
            logical_readout = (logical_readout + readout_values[id]) % 2

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

    def logical_x_plus_reset(self):
        self.circ.reset(self.qregisters["data"])
        self.circ.h(self.qregisters["data"][1::2])  # H|0> = |+>

    def logical_z_plus_reset(self):
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
        # add H to measure along X for odd data qubits
        self.circ.h(self.qregisters["data"][0::2])

        self.circ.measure(self.qregisters["data"], self.cregisters["lattice_readout"])
        self.circ.barrier()

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
        Returns:
            logical_readout (int): 
                logical readout value
            syndromes (Dict[str, List[TQubit]]]):
                key: syndrome type
                value: (time, row, col) of parsed syndrome hits (changes between consecutive rounds)
        """
        chunks = readout_string.split(" ")

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

        mask_B = "1" * self.params["num_syn"]
        mask_A = mask_B + "0" * self.params["num_syn"]
        A_syndromes = [
            (x & int(mask_A, base=2)) >> self.params["num_syn"] for x in xor_syndromes
        ]
        B_syndromes = [x & int(mask_B, base=2) for x in xor_syndromes]

        A = []
        for T, syndrome in enumerate(A_syndromes):
            for loc in range(self.params["num_syn"]):
                if syndrome & 1 << loc:
                    A.append((float(T), -0.5 + loc, 0.5 + loc % 2))

        B = []
        for T, syndrome in enumerate(B_syndromes):
            for loc in range(self.params["num_syn"]):
                if syndrome & 1 << loc:
                    B.append((float(T), 0.5 + loc // 2, 0.5 + loc % 2 * 2 - loc // 2))

        return (
            logical_readout,
            {"A": A, "B": B,},
        )


class XZZXQubit(TopologicalQubit[TQubit]):
    """
    A single logical surface code qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend TopologicalQubit which extends QuantumCircuit.
    """

    def __init__(
        self,
        params: Dict[str, int],
        name: str = "tqubit",
        circ: Optional[QuantumCircuit] = None,
    ) -> None:
        """
        Initializes a new QuantumCircuit for this logical qubit and calculates
        the underlying surface code lattice ordering.
        
        Args:
            d (int): Number of physical "data" qubits. Only odd d is possible!
        """
        # == None is necessary, as "not circ" is true for circ=QuantumCircuit()
        circ = QuantumCircuit() if circ == None else circ

        super().__init__(circ, name)
        self.lattice = _XZZXLattice(params, name, circ)

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
            self.lattice.qregisters["mb"],
            syndrome_readouts[0 : self.lattice.params["num_syn"]],
        )
        self.circ.measure(
            self.lattice.qregisters["ma"],
            syndrome_readouts[
                self.lattice.params["num_syn"] : self.lattice.params["num_syn"] * 2
            ],
        )
        self.circ.reset(self.lattice.qregisters["mb"])
        self.circ.reset(self.lattice.qregisters["ma"])
        self.circ.barrier()

    def identity(self) -> None:
        """
        Inserts an identity on the data and syndrome qubits. This is a hack to
        create an isolated error model.
        """
        [self.circ.id(register) for register in self.lattice.qregisters.values()]
        self.circ.barrier()

    def identity_data(self) -> None:
        """
        Inserts an identity on the data qubits only. This is a hack to create an
        isolated error model.
        """
        self.circ.id(self.lattice.qregisters["data"])
        self.circ.barrier()

    def logical_x_plus_reset(self):
        """
        Initialize/reset to a logical |x+> state.
        """
        self.lattice.logical_x_plus_reset()

    def logical_z_plus_reset(self):
        """
        Initialize/reset to a logical |z+> state.
        """
        self.lattice.logical_z_plus_reset()

    def logical_x(self) -> None:
        """
        Logical X operator on the topological qubit.
        Defined as the left-most column on the A Syndrome Graph.
        """
        self.lattice.logical_x()

    def logical_z(self) -> None:
        """
        Logical Z operator on the topological qubit.
        Defined as the top-most row on the B Syndrome Graph.
        """
        self.lattice.logical_z()

    def readout_x(self) -> None:
        """
        Convenience method to read-out the logical-X projection.
        Defined as the left-most column.
        """
        self.lattice.readout_x()

    def readout_z(self) -> None:
        """
        Convenience method to read-out the logical-Z projection.
        Defined as the top-most row.
        """
        self.lattice.readout_z()

    def lattice_readout_x(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of A stabilizer measurments,
        as well as a logical X readout.
        """

        self.lattice.lattice_readout_x()

    def lattice_readout_z(self) -> None:
        """
        Readout all data qubits that constitute the lattice.
        This readout can be used to extract a final round of B stabilizer measurments,
        as well as a logical Z readout.
        """

        self.lattice.lattice_readout_z()

    def parse_readout(
        self, readout_string: str, readout_type: Optional[str] = None
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
        return self.lattice.parse_readout(readout_string, readout_type)
