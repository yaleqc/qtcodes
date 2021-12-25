"""
Repetition Code Encoder Classes
"""
from numbers import Number
from typing import Dict, Tuple, Optional
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
        params = self._validate_params(params)
        super().__init__(params, name, circ)

    @staticmethod
    def _validate_params(params: Optional[Dict[str, int]]) -> Dict[str, int]:
        params = params if params else {}

        if "phase-flip-protected" not in params:
            params["phase-flip-protected"] = False

        # The repetition code is a 1D lattice
        if "d" not in params:
            if params["phase-flip-protected"]:
                params["d"] = (1, 3)
            else:
                params["d"] = (3, 1)
        elif isinstance(params["d"], Number):
            d = int(params["d"])
            if params["phase-flip-protected"]:
                params["d"] = (1, d)
            else:
                params["d"] = (d, 1)
        elif isinstance(params["d"], Tuple):
            if not params["phase-flip-protected"] and params["d"][constants.DW] != 1:
                raise LatticeError(
                    "Bit-flip protected repetition qubits can only have width 1 in parameter d: e.g. (3,1). If you "
                    + "intend to create a phase-flip protected repetition qubit (e.g. with d=(1,3)), then please set the "
                    + "phase-flip-protected flag parameter to True."
                )
            if params["phase-flip-protected"] and params["d"][constants.DH] != 1:
                raise LatticeError(
                    "Phase-flip protected repetition qubits can only have height 1 in parameter d: e.g. (1,3). If you "
                    + "intend to create a bit-flip protected repetition qubit (e.g. with d=(1,3)), then please set the "
                    + "phase-flip-protected flag parameter to False."
                )
        else:
            raise LatticeError(
                "Please specify either an integer value (e.g. 3) or tuple (e.g. (3, 1)) for parameter d. If you intend "
                + "to create a phase-flip protected repetition qubit (e.g. with d=(1,3)), then please set the "
                + "phase-flip-protected flag parameter to True."
            )

        return params
