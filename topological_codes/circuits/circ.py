from topological_codes import XXZZQubit, RepetitionQubit, TopologicalQubit, XZZXQubit
from qiskit import QuantumCircuit
from typing import Union, Dict


class TopologicalRegister:
    """
    A blueprint for a TopologicalRegister that stores topological qubit(s)
    """
    def __init__(
        self,
        num_tqubits: int,
        circ: QuantumCircuit = None,
        type: str = "XXZZ",
        params: Dict[str, int] = {"d": 3},
        name: str = "treg",
    ):
    
        """
        Args:
            num_tqubits: int
                The number of topological qubits to be stored in the register.
            circ: QuantumCircuit
                QuantumCircuit on top of which the topological qubit is built. 
                This is often shared amongst multiple TQubits.
                If none is provided, then a new QuantumCircuit is initialized and stored.
             type: str
                Specifies the type of TopologicalQubit
                Must be a value of the blueprint dictionary
             params (Dict[str,int]): 
                 Contains params such as d, where d is the number of physical "data" qubits lining a row or column of the lattice. 
                 Only odd d is possible!    
             name: str
                Useful when combining multiple TopologicalQubits together. Prepended to all registers.
            
        """
        self.name = name

        # == None is necessary, as "not circ" is true for circ=QuantumCircuit()
        self.circ = QuantumCircuit() if circ == None else circ

        self.tqubits = []

        
        blueprint = {"XXZZ": XXZZQubit, "Repetition": RepetitionQubit, "XZZX": XZZXQubit}
        if type not in blueprint:
            raise ValueError(
                "Please choose a Topological Qubit type from: "
                + str(list(blueprint.keys()))
            )
        self.tqubit_type = blueprint[type]

        for i in range(num_tqubits):
            self.tqubits.append(
                self.tqubit_type(
                    params=params, name=self.name + "_" + str(i), circ=self.circ
                )
            )
    
    def __getitem__(self, key: int):
        """
        Allows us to return the nth element of TopologicalRegister as a list. 
        """
        return self.tqubits[key]


class TopologicalCircuit: 
    """
    A blueprint for a TopologicalCircuit. 
    Shares the same QuantumCircuit object created in TopologicalRegister.
    """
    def __init__(self, treg: TopologicalRegister):
        self.treg = treg
        self.circ = treg.circ 

    def _get_index(self, tqubit: Union[TopologicalQubit, int]):
        if type(tqubit) == int:
            tqubit = self.treg[tqubit]
        assert isinstance(tqubit, TopologicalQubit)
        return tqubit

    def logical_x_plus_reset(self, tqubit: Union[TopologicalQubit, int]):
        """
        Initialize/reset to a logical |x+> state.
        """
        tqubit = self._get_index(tqubit)
        tqubit.logical_x_plus_reset()

    def logical_z_plus_reset(self, tqubit: Union[TopologicalQubit, int]):
        """
        Initialize/reset to a logical |z+> state.
        """
        tqubit = self._get_index(tqubit)
        tqubit.logical_z_plus_reset()
     
    def logical_x(self, tqubit: Union[TopologicalQubit, int]):
        """
        Applies a logical X operator
        """
        tqubit = self._get_index(tqubit)
        tqubit.logical_x()
        
    def logical_z(self, tqubit: Union[TopologicalQubit, int]):
        """
        Applies a logical Z operator
        """
        tqubit = self._get_index(tqubit)
        tqubit.logical_z()
        
    def stabilize(self, tqubit: Union[TopologicalQubit, int]):
        """
        Runs a single round of stabilization (entangle and measure)
        """ 
        tqubit = self._get_index(tqubit)
        tqubit.stabilize()

    def identity(self, tqubit: Union[TopologicalQubit, int]):
        """
        Inserts an identity on the data and syndrome qubits. This is a hack to
        create an isolated error model.
        """
        tqubit = self._get_index(tqubit)
        tqubit.identity()
        
    def identity_data(self, tqubit: Union[TopologicalQubit, int]):
        """
        Inserts an identity on the data qubits only. This is a hack to create an
        isolated error model.
        """
        tqubit = self._get_index(tqubit)
        tqubit.identity_data()

    def hadamard_reset(self, tqubit: Union[TopologicalQubit, int]):
        """
        A hack to initialize a + and - logical qubit for now...
        """
        tqubit = self._get_index(tqubit)
        tqubit.hadamard_reset()
        
    def readout_z(self, tqubit: Union[TopologicalQubit, int]):
        """
        Convenience method to read-out the logical-Z projection.
        """
        tqubit = self._get_index(tqubit)
        tqubit.readout_z()
        
    def readout_x(self, tqubit: Union[TopologicalQubit, int]):
        """
        Convenience method to read-out the logical-X projection.
        """
        tqubit = self._get_index(tqubit)
        tqubit.readout_x()
    
    def draw(self, **kwargs):
        return self.circ.draw(**kwargs)

    def __str__(self):
        return self.circ.__str__()
    
    def draw_simple():
        #TODO: Gate that spans multiple qubits when we draw our circuit
        pass
    


