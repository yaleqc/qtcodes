from abc import abstractmethod, ABCMeta
from typing import Tuple, List, Dict, Optional, TypeVar, Generic, Union
import networkx as nx

TQubit = TypeVar("TQubit")


class TopologicalGraphDecoder(Generic[TQubit], metaclass=ABCMeta):
    """
    Abstract class for topological code MWPM decoders to implement.
    """

    @abstractmethod
    def __init__(self, code_params: Dict) -> None:
        self.code_params = code_params
        self.S: Dict[str, nx.Graph] = {}  # syndrome graphs

    @abstractmethod
    def _make_syndrome_graph(self) -> None:
        pass

    @abstractmethod
    def correct_readout(
        self,
        syndromes: Union[str, Dict[str, List[TQubit]]],
        logical_qubit_value: Optional[int] = None,
        logical_readout_type: str = "Z",
    ) -> int:
        pass

    @abstractmethod
    def _convert_string_to_nodes(
        self, readout_string: str
    ) -> Tuple[int, Dict[str, List[TQubit]]]:
        pass

    @abstractmethod
    def _make_error_graph(
        self, nodes: List[TQubit], syndrome_graph_key: str, err_prob: Optional[int]
    ) -> nx.Graph:
        pass

    @abstractmethod
    def _make_matching_graph(
        self, error_graph: nx.Graph, syndrome_graph_key: str
    ) -> nx.Graph:
        pass

    @abstractmethod
    def _run_mwpm(
        self, matching_graph: nx.Graph, syndrome_graph_key: str,
    ) -> List[Tuple[TQubit, TQubit]]:
        pass
