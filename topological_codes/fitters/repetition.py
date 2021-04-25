# -*- coding: utf-8 -*-
"""
Graph decoder for surface codes
"""
import copy
import math
from itertools import combinations, product
from collections import defaultdict

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from typing import Tuple, List, Dict, Optional, TypeVar, Union, cast
from .base import TopologicalGraphDecoder
from .visualization import VisualizationMixin
from .mwpm import MWPMDecodingMixin


from qiskit import QuantumCircuit, execute

try:
    from qiskit import Aer

    HAS_AER = True
except ImportError:
    from qiskit import BasicAer

    HAS_AER = False


TQubit = Tuple[float, float, float]  # (time,row,column) ==> (t,i,j)
TQubitLoc = Tuple[float, float]  # (row,column) ==> (i,j)


class RepetitionGraphDecoderBase(TopologicalGraphDecoder[TQubit]):
    def __init__(self, code_params: Dict) -> None:
        self.code_params = code_params
        self.S: Dict[str, nx.Graph] = {}  # syndrome graphs

    def _make_syndrome_graph(self) -> None:
        pass

    def correct_readout(
        self,
        syndromes: Union[str, Dict[str, List[TQubit]]],
        logical_qubit_value: Optional[int] = None,
        logical_readout_type: str = "Z",
    ) -> int:
        pass

    def _string2nodes(self, readout_string: str) -> Tuple[int, Dict[str, List[TQubit]]]:
        pass

    def _make_error_graph(
        self, nodes: List[TQubit], syndrome_graph_key: str, err_prob: Optional[int]
    ) -> nx.Graph:
        pass

    def _run_mwpm(
        self, matching_graph: nx.Graph, syndrome_graph_key: str,
    ) -> List[Tuple[TQubit, TQubit]]:
        pass
