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

        self.look_up_table: Dict[str, Dict[int, Dict[str, int]]] = {}
        self.look_up_table_shots: Dict[str, Dict[int, int]] = {}
        self.reset_look_up_table()

    @abstractmethod
    def _make_syndrome_graph(self) -> None:
        """
        Construct syndrome graphs used for MWPM decoding.
        """

    @abstractmethod
    def correct_readout(
        self,
        syndromes: Union[str, Dict[str, List[TQubit]]],
        logical_readout_type: str,
        logical_qubit_value: Optional[int] = None,
        err_prob: Optional[float] = None,
    ) -> int:
        """
        Args:
            syndromes (Union[str, Dict[str, List[TQubit]]]): either...
                (Dict[str, List[TQubit]]]):
                    key (str): syndrome type
                    value (TQubit): syndrome node hits (changes between consecutive rounds)
                (str): readout string
            logical_qubit_value (Optional[int]): measured logical qubit value
            logical_readout_type (str): logical readout type (e.g. "X" or "Z")
            err_prob (Optional[float]): Probability of IID data qubit X/Z flip. Defaults to None.
        Returns:
            logical_qubit_value (int): The most probable encoded value of the logical qubit.
        Additional Information:
            This method can be used to benchmark logical error rates, as well as perform fault tolerant readout.
        """

    @abstractmethod
    def parse_readout(self, readout_string: str) -> Tuple[int, Dict[str, List[TQubit]]]:
        """
        Converte between readout_string to logical_readout and syndrome nodes.

        Args:
            readout_string (str): readout string from quantum circuit

        Returns:
            logical_readout (int): logical readout value from readout
            syndromes (Dict):
                key (str): syndrome type key
                val (List[TQubit]): list of syndrome nodes hit
        """

    @abstractmethod
    def _make_error_graph(
        self,
        nodes: List[TQubit],
        syndrome_graph_key: str,
        err_prob: Optional[float] = None,
    ) -> rx.PyGraph:
        """Creates error syndrome subgraph from list of syndrome nodes. The output of
        this function is a graph that's ready for minimum weight perfect matching (MWPM).

        If err_prob is specified, we adjust the shortest distance between syndrome
        nodes by the degeneracy of the error path.

        Args:
            nodes ([TQubit,]): List of changes of syndrome nodes in time.
            syndrome_graph_key (char): Which syndrome subgraph these nodes are from.
            err_prob (Optional[float]): Probability of IID data qubit X/Z flip. Defaults to None.

        Returns:
            error_graph (rx.PyGraph): Nodes are syndromes, edges are proxy for error probabilities
        """

    @abstractmethod
    def _run_mwpm(self, matching_graph: rx.PyGraph) -> List[Tuple[TQubit, TQubit]]:
        """
        Return matches of minimum weight perfect matching (MWPM) on matching_graph.

        Args:
            matching_graph (rx.PyGraph): Graph to run MWPM.
            floats (bool):
                whether or not the graph contains float edge weights

        Returns:
            [(TQubit, TQubit),]: List of matchings found from MWPM
        """

    def reset_look_up_table(self) -> None:
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
        self.look_up_table = {
            "X": {0: {}, 1: {}},
            "Z": {0: {}, 1: {}},
        }
        self.look_up_table_shots = {
            "X": {0: 0, 1: 0},
            "Z": {0: 0, 1: 0},
        }

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

    def _run_look_up(
        self,
        syndromes: str,
        logical_readout_type: str,
        p_i: Optional[List[float]] = None,
    ) -> List[float]:
        """
        Calculates conditional prob p(|i> | syndromes) by look up table decoding.

        Args:
            syndromes (str):
                Readout of the form "0 00000000 00000000" (logical_readout syndrome_1 syndrome_0)
                or of the form "000000000 00000000 00000000" (lattice_readout syndrome_1 syndrome_0)
            logical_readout_type (str): "X" or "Z"
            p_i (Optional[List[float]]): p_i[i] is the posterior probability that the logical qubit is |i>

        Returns:
            p_i_given_s (List[float]):
                p(i|s) conditional probability that given syndromes s, the qubit is in |i>
        """

        p_s_given_i = [
            0.0,
            0.0,
        ]  # p_s_given_i[i] is the prob that syndromes s was produced by |i>, i.e. p(s|i)
        for log in [0, 1]:
            if syndromes in self.look_up_table[logical_readout_type][log]:
                p_s_given_i[log] = (
                    self.look_up_table[logical_readout_type][log][syndromes]
                    / self.look_up_table_shots[logical_readout_type][log]
                )

        # posterior prob that the qubit in |i>, i.e. p(i)
        p_i = [1.0 / 2, 1.0 / 2] if not p_i else p_i

        # p(s) = sum_i p(s|i)p(i)
        p_s = p_s_given_i[0] * p_i[0] + p_s_given_i[1] * p_i[1]

        # p(i|s) = p(s|i)p(i)/p(s), Bayes Rule
        p_i_given_s = [p_s_given_i[i] * p_i[i] / p_s for i in [0, 1]]
        return p_i_given_s

    def correct_readout_look_up_table(
        self, syndromes: str, logical_readout_type: str
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
        p_i_given_s = self._run_look_up(syndromes, logical_readout_type)
        return p_i_given_s.index(max(p_i_given_s))

