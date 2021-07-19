"""
Quantum Topological Codes (qtcodes) Package
"""
import os

from .circuits import (
    RepetitionQubit,
    XXZZQubit,
    XZZXQubit,
    TopologicalRegister,
    TopologicalCircuit,
)

from .fitters import RotatedDecoder, RepetitionDecoder

from .tools import TopologicalBenchmark, TopologicalAnalysis, TopologicalBatchAnalysis

with open(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "VERSION.txt")), "r"
) as _ver_file:
    __version__ = _ver_file.read().rstrip()

__author__ = "Shantanu Jha"
__credits__ = "Qiskit Topological Codes Team"
