"""
Quantum Topological Codes (qtcodes) Package
"""
import os

from .circuits import *
from .fitters import *
from .tools import *

with open(
    os.path.abspath(os.path.join(os.path.dirname(__file__), "VERSION.txt")), "r"
) as _ver_file:
    __version__ = _ver_file.read().rstrip()

__author__ = "Shantanu Jha"
__credits__ = "Qiskit Topological Codes Team"
