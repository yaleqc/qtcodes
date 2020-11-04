# -*- coding: utf-8 -*-
"""
benchmarking class for surface codes
"""
import numpy as np
import matplotlib.pyplot as plt
from .fitters import GraphDecoder


class SurfaceCodeBenchmarkingTool:
    def __init__(self, decoder):
        self.decoder = decoder

    def logical_error_rate(self, readout_strings, correct_logical_value):
        """
        Args:
            decoder: a GraphDecoder object corresponding to the surface code used to produce the readout_strings

            readout_strings: a dictionary of readout strings along with counts
            e.g. {"1 00000000 00000000":48, "1 00100000 00100000":12, ...} in the case of d=3 and T=2

            correct_logical_value: integer (0/1) depicting original encoded logical value 

        Returns:
            error_rate: float = (number of unsuccessful logical value predictions) / (total number of predictions )
        """
        total_count = 0
        total_errors = 0
        for readout, count in readout_strings.items():
            total_count += count
            predicted_logical_value = self.decoder.correct_readout(readout)
            if predicted_logical_value != correct_logical_value:
                total_errors += count

        return total_errors / total_count


"""
    def simulate_readout(decoder, correct_logical_value):
        return readout_strings
"""
