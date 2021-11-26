# -*- coding: utf-8 -*-
"""
Graph decoder for rep code
"""
from numbers import Number

from qtcodes.circuits.repetition import RepetitionQubit
from qtcodes.fitters.rotated_surface import RotatedDecoder


class RepetitionDecoder(RotatedDecoder):
    """
    Class to construct the graph corresponding to the possible syndromes
    of a quantum error correction Repetition code, and then run suitable decoders.
    """

    encoder_type = RepetitionQubit
    syndrome_graph_keys = ["Z"]

    def _params_validation(self):
        # The repetition code is a 1D lattice
        if "d" not in self.params:
            self.params["d"] = (3, 1)
        elif isinstance(self.params["d"], Number):
            d = int(self.params["d"])
            self.params["d"] = (d, 1)
        else:
            self.params["d"] = (self.params["d"][0], 1)

        super()._params_validation()
