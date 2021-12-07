"""
Repetition Code Encoder Classes
"""
from numbers import Number
from typing import Dict, List, Tuple, Optional
from qiskit import QuantumCircuit

from qtcodes.circuits.xxzz import XXZZQubit
from qtcodes.circuits.base import LatticeError
from qtcodes.common import constants

TQubit = Tuple[float, float, float]


class RepetitionQubit(XXZZQubit):
    """
    A single logical repetition code qubit. At the physical level, this wraps a
    circuit, so we chose to subclass and extend TopologicalQubit which extends QuantumCircuit.
    """

    def __init__(
        self,
        params: Optional[Dict[str, int]] = None,
        name: str = "tq",
        circ: Optional[QuantumCircuit] = None,
    ) -> None:
        params = params if params else {}

        # The repetition code is a 1D lattice
        if "d" not in params:
            params["d"] = (3, 1)
        elif isinstance(params["d"], Number):
            d = int(params["d"])
            params["d"] = (d, 1)
        elif isinstance(params["d"], Tuple):
            if params["d"][constants.DW] != 1:
                raise LatticeError(
                    "Repetition qubits can only have width 1 in parameter d: e.g. (3,1)."
                )
        else:
            raise LatticeError(
                "Please provide a valid height in parameter d: e.g. 3 or (3,1)."
            )

        super().__init__(params, name, circ)
