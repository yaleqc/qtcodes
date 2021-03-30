from qiskit import QuantumRegister, QuantumCircuit


class TopologicalEntangler:
    def __init__(self, circ, stabilizers, qubit_indices_lists):
        self.circ = circ
        self.stabilizers = stabilizers
        self.qubit_indices_lists = qubit_indices_lists

    def entangle(self):
        for i in range(len(self.stabilizers)):
            stabilizer = self.stabilizers[i](self.circ, self.qubit_indices_lists[i])
            stabilizer.entangle()


class _Parity:
    def __init__(self, circ, qubit_indices):
        self.circ = circ
        self.qubit_indices = qubit_indices

    def entangle(self):
        link_qubit = self.qubit_indices[0]
        code_qubit_0 = self.qubit_indices[1]
        code_qubit_1 = self.qubit_indices[2]
        self.circ.cx(code_qubit_0, link_qubit)
        self.circ.cx(code_qubit_1, link_qubit)


class RepetitionCode:
    def __init__(self, d):
        self.d = d
        self.circ = QuantumCircuit(2 * d - 1, 2 * d - 1)

    def stabilize(self):
        stabilizers = [_Parity, _Parity]
        qubit_indices_lists = []
        for i in range(self.d):
            qubit_indices_lists.append([i, self.d - 1 + i, self.d + i])
        entangler = TopologicalEntangler(self.circ, stabilizers, qubit_indices_lists)
        entangler.entangle()
        self.circ.barrier()
        self.circ.measure(list(range(self.d)), list(range(self.d)))

