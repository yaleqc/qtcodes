"""
Base Topological Decoder Classes
"""
from abc import abstractmethod, ABCMeta
from typing import Tuple, List, Dict, Optional, TypeVar, Generic, Union, Type
import retworkx as rx

TQubit = TypeVar("TQubit")


class TopologicalDecoder(Generic[TQubit], metaclass=ABCMeta):
    """
    Abstract class for topological code MWPM decoders to implement.
    """

    @property
    @abstractmethod
    def encoder_type(self) -> Type:
        """
        TopologicalQubit
        """

    @abstractmethod
    def __init__(self, params: Dict) -> None:
        self.S: Dict[str, rx.PyGraph] = {}
        self.node_map: Dict[str, Dict[TQubit, int]] = {}
        self.params = params

        table, shots = self.empty_look_up_table()
        self.look_up_table: Dict[str, Dict[int, Dict[str, int]]] = table
        self.look_up_table_shots: Dict[str, Dict[int, int]] = shots

    @abstractmethod
    def _make_syndrome_graph(self) -> None:
        pass

    @abstractmethod
    def correct_readout(
        self,
        syndromes: Union[str, Dict[str, List[TQubit]]],
        logical_readout_type: str,
        logical_qubit_value: Optional[int] = None,
        err_prob: Optional[float] = None,
    ) -> int:
        pass

    @abstractmethod
    def _string2nodes(self, readout_string: str) -> Tuple[int, Dict[str, List[TQubit]]]:
        pass

    @abstractmethod
    def _make_error_graph(
        self,
        nodes: List[TQubit],
        syndrome_graph_key: str,
        err_prob: Optional[float] = None,
    ) -> rx.PyGraph:
        pass

    @abstractmethod
    def _run_mwpm(self, matching_graph: rx.PyGraph) -> List[Tuple[TQubit, TQubit]]:
        pass

    def empty_look_up_table(
        self,
    ) -> Tuple[Dict[str, Dict[int, Dict[str, int]]], Dict[str, Dict[int, int]]]:
        """
        Empty look up table used to reset the self.look_up_table

        Returns:
            (Dict[str, Dict[int, Dict[str, int]]]):
                key (str): logical_readout_type
                val: key (int): logical_qubit_value
                     val: key (str): readout_string
                          val (int): counts
            (Dict[str, Dict[int, int]]):
                key (str): logical_readout_type
                val: key (int): logical_qubit_value
                     val (int): total_counts
        """

        return (
            {"X": {0: {}, 1: {}}, "Z": {0: {}, 1: {}},},
            {"X": {0: 0, 1: 0}, "Z": {0: 0, 1: 0},},
        )

    def set_look_up_table(
        self,
        logical_readout_type: str,
        logical_qubit_value: int,
        results: Dict[str, int],
    ):
        """
        Used to set look up table.

        Args:
            logical_readout_type (str): "X" or "Z"
            logical_qubit_value (int): 0 or 1
            results:
                key (str): readout_string
                val (int): counts
        """
        self.look_up_table[logical_readout_type][logical_qubit_value] = results
        self.look_up_table_shots[logical_readout_type][logical_qubit_value] = sum(
            results.values()
        )

    def correct_readout_look_up_table(
        self, syndromes: str, logical_readout_type: str,
    ) -> int:
        """
        Calculates most likely logical_qubit_value by look up table decoding.

        Args:
            syndromes (str):
                Readout of the form "0 00000000 00000000" (logical_readout syndrome_1 syndrome_0)
                or of the form "000000000 00000000 00000000" (lattice_readout syndrome_1 syndrome_0)
            logical_readout_type (str): "X" or "Z"

        Returns:
            logical_qubit_value (int):
                most probable original logical qubit value from look up table decoding
        """
        p = [0.0, 0.0]
        for log in [0, 1]:
            if syndromes in self.look_up_table[logical_readout_type][log]:
                p[log] = (
                    self.look_up_table[logical_readout_type][log][syndromes]
                    / self.look_up_table_shots[logical_readout_type][log]
                )
        return p.index(max(p))

