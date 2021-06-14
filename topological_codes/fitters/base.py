"""
Base Topological Decoder Classes
"""
from abc import abstractmethod, ABCMeta
from typing import Tuple, List, Dict, Optional, TypeVar, Generic, Union, Type
import retworkx as rx

TQubit = TypeVar("TQubit")


class TopologicalGraphDecoder(Generic[TQubit], metaclass=ABCMeta):
    """
    Abstract class for topological code MWPM decoders to implement.
    """

    S: Dict[str, rx.PyGraph] = {}
    node_map: Dict[str, Dict[TQubit, int]] = {}

    @property
    @abstractmethod
    def encoder_type(self) -> Type:
        """
        TopologicalQubit
        """

    @abstractmethod
    def __init__(self, params: Dict) -> None:
        self.params = params

    @abstractmethod
    def _make_syndrome_graph(self) -> None:
        pass

    @abstractmethod
    def correct_readout(
        self,
        syndromes: Union[str, Dict[str, List[TQubit]]],
        logical_qubit_value: Optional[int] = None,
        logical_readout_type: str = "Z",
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
