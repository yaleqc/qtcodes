# -*- coding: utf-8 -*-
"""
Graph decoder for rep code
"""
from numbers import Number
from typing import Tuple

from qtcodes.circuits.repetition import RepetitionQubit
from qtcodes.fitters.rotated_surface import RotatedDecoder
from qtcodes.circuits.base import LatticeError
from qtcodes.common import constants


class RepetitionDecoder(RotatedDecoder):
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction Repetition code, and then run suitable decoders.
    """

    encoder_type = RepetitionQubit
    syndrome_graph_keys = ["Z"]

    def _params_validation(self):
        if "phase-flip-protected" not in self.params:
            self.params["phase-flip-protected"] = False

        if self.params["phase-flip-protected"]:
            self.syndrome_graph_keys = ["X"]

        # The repetition code is a 1D lattice
        if "d" not in self.params:
            if self.params["phase-flip-protected"]:
                self.params["d"] = (1, 3)
            else:
                self.params["d"] = (3, 1)
        elif isinstance(self.params["d"], Number):
            d = int(self.params["d"])
            if self.params["phase-flip-protected"]:
                self.params["d"] = (1, d)
            else:
                self.params["d"] = (d, 1)
        elif isinstance(self.params["d"], Tuple):
            if (
                not self.params["phase-flip-protected"]
                and self.params["d"][constants.DW] != 1
            ):
                raise LatticeError(
                    "Repetition qubits can only have width 1 in parameter d: e.g. (3,1)."
                )
            if (
                self.params["phase-flip-protected"]
                and self.params["d"][constants.DH] != 1
            ):
                raise LatticeError(
                    "Phase-flip protected repetition qubits can only have height 1 in parameter d: e.g. (1,3)."
                )
        else:
            if self.params["phase-flip-protected"]:
                raise LatticeError(
                    "Please provide a valid width in parameter d: e.g. 3 or (1,3)."
                )
            else:
                raise LatticeError(
                    "Please provide a valid height in parameter d: e.g. 3 or (3,1)."
                )

        super()._params_validation()
