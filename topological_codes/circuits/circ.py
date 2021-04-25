from topological_codes import XXZZQubit, RepetitionQubit, TopologicalQubit
from qiskit import QuantumCircuit
from typing import Union


class TopologicalRegister:
    def __init__(
        self, num_tqubits: int, circ=None, type="XXZZ", params={"d": 3}, name="treg"
    ):
        self.name = name

        # == None is necessary, as "not circ" is true for circ=QuantumCircuit()
        self.circ = QuantumCircuit() if circ == None else circ

        self.tqubits = []

        blueprint = {"XXZZ": XXZZQubit, "Repetition": RepetitionQubit}
        if type not in blueprint:
            raise ValueError(
                "Please choose a Topological Qubit type from: "
                + str(list(blueprint.keys()))
            )
        self.tqubit_type = blueprint[type]

        for i in range(num_tqubits):
            self.tqubits.append(
                self.tqubit_type(params=params, name=self.name + str(i), circ=self.circ)
            )

    def __getitem__(self, key: int):
        return self.tqubits[key]


class TopologicalCircuit:
    def __init__(self, treg: TopologicalRegister):
        self.treg = treg
        self.circ = treg.circ

    def get_index(self, tqubit: Union[TopologicalQubit, int]):
        if type(tqubit) == int:
            tqubit = self.treg[tqubit]
        assert isinstance(tqubit, TopologicalQubit)
        return tqubit

    def x(self, tqubit: Union[TopologicalQubit, int]):
        tqubit = self.get_index(tqubit)
        tqubit.logical_x()

    def stabilize(self, tqubit: Union[TopologicalQubit, int]):
        tqubit = self.get_index(tqubit)
        tqubit.stabilize()

    def draw(self, **kwargs):
        return self.circ.draw(**kwargs)

    def __str__(self):
        return self.circ.__str__()
