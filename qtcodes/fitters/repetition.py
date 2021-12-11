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
        self.params = RepetitionQubit._validate_params(self.params)

        if self.params["phase-flip-protected"]:
            self.syndrome_graph_keys = ["X"]

        super()._params_validation()
