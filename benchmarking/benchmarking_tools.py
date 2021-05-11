# -*- coding: utf-8 -*-
"""
benchmarking class for surface codes
"""
import sys
import os

sys.path.insert(0, ".." + os.sep)

import numpy as np
import matplotlib.pyplot as plt
from topological_codes import RotatedGraphDecoder

from qiskit import QuantumCircuit, execute, QuantumRegister, ClassicalRegister, Aer
from tqdm import tqdm


plt.rcParams.update({"font.size": 14, "pdf.fonttype": 42, "ps.fonttype": 42})

from multiprocessing import Pool


class SurfaceCodeBenchmarkingTool:
    def __init__(
        self,
        decoder=None,
        readout_circuit=None,
        noise_model_func=None,
        filename=None,
        correct_logical_value=None,
    ):
        self.decoder = decoder
        if self.decoder is not None:
            self.d = decoder.params["d"]
            self.T = decoder.params["T"]
        self.filename = (
            "surface_code_d_{}_T_{}.npz".format(int(self.d), int(self.T))
            if filename is None
            else filename
        )
        self.readout_circuit = readout_circuit
        self.correct_logical_value = correct_logical_value
        self.noise_model_func = noise_model_func
        self.benchmark_data = {"noise": [], "logical_error_rate": []}

    def logical_error_rate(self, readout_strings, correct_logical_value, err_prob=None):
        """
        Args:
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
            predicted_logical_value = self.decoder.correct_readout(
                readout, err_prob=err_prob
            )
            if predicted_logical_value != correct_logical_value:
                total_errors += count

        return total_errors / total_count

    def simulate_readout(
        self,
        correct_logical_value=0,
        noise_values=[
            5e-5,
            1e-4,
            2e-4,
            5e-4,
            1e-3,
            2e-3,
            4e-3,
            5e-3,
            6e-3,
            7e-3,
            8e-3,
            9e-3,
            1e-2,
            2e-2,
        ],
        save_data=True,
        shots=1024 * 2,
    ):
        self.benchmark_data["noise"] = []
        self.benchmark_data["logical_error_rate"] = []

        noise_values = sorted(
            noise_values, reverse=True
        )  # higher noise readout is slower to decode, gives more accurate tqdm estimate
        for noise_value in tqdm(noise_values):
            results = (
                execute(
                    self.readout_circuit,
                    Aer.get_backend("qasm_simulator"),
                    noise_model=self.noise_model_func(noise_value),
                    shots=shots,
                )
                .result()
                .get_counts()
            )
            logical_error_rate_value = self.logical_error_rate(
                results, correct_logical_value
            )
            self.benchmark_data["noise"] = [noise_value] + self.benchmark_data["noise"]
            self.benchmark_data["logical_error_rate"] = [
                logical_error_rate_value
            ] + self.benchmark_data["logical_error_rate"]
            if save_data:
                self.save_data()

        return self.benchmark_data

    def simulate_readout_mp(
        self,
        correct_logical_value=0,
        noise_values=[
            5e-5,
            1e-4,
            2e-4,
            5e-4,
            1e-3,
            2e-3,
            4e-3,
            5e-3,
            6e-3,
            7e-3,
            8e-3,
            9e-3,
            1e-2,
            2e-2,
        ],
        save_data=True,
    ):
        noise_values = sorted(noise_values)
        self.correct_logical_value = correct_logical_value
        p = Pool(4)  # TODO change on HPC

        logical_error_rates = p.map(self.simulate_readout_single, noise_values)

        self.benchmark_data["noise"] = noise_values
        self.benchmark_data["logical_error_rate"] = logical_error_rates
        self.save_data()
        return self.benchmark_data

    def simulate_readout_single(self, noise_value):
        results = (
            execute(
                self.readout_circuit,
                Aer.get_backend("qasm_simulator"),
                noise_model=self.noise_model_func(noise_value),
                shots=100000,
            )
            .result()
            .get_counts()
        )
        logical_error_rate_value = self.logical_error_rate(
            results, self.correct_logical_value, err_prob=noise_value
        )
        print("Done simulating noise: " + str(noise_value))
        self.benchmark_data["noise"].append(noise_value)
        self.benchmark_data["logical_error_rate"].append(logical_error_rate_value)
        self.save_data()
        return logical_error_rate_value

    def plot_benchmark_data(
        self, fig=None, ax=None, log=True, per_round=False, **kwargs
    ):
        if fig is None:
            fig = plt.figure(figsize=(3.5, 2.5), dpi=200)
        if ax is None:
            ax = fig.subplots()
        plt.plot(
            self.benchmark_data["noise"],
            np.array(self.benchmark_data["logical_error_rate"])
            / (self.T if per_round else 1.0),
            **kwargs
        )
        if log:
            plt.yscale("log")
            plt.xscale("log")

    def save_data(self, filename=None):
        filename = self.filename if filename is None else filename
        np.savez(
            filename,
            d=self.decoder.params["d"],
            T=self.decoder.params["T"],
            noise=self.benchmark_data["noise"],
            logical_error_rate=self.benchmark_data["logical_error_rate"],
        )

    def load_data(self, filename=None):
        filename = self.filename if filename is None else filename
        data = np.load(filename)
        self.d = int(data["d"])
        self.T = int(data["T"])
        self.decoder = RotatedGraphDecoder({"d": self.d, "T": self.T})

        # self.readout_circuit =
        self.benchmark_data["noise"] = data["noise"]
        self.benchmark_data["logical_error_rate"] = data["logical_error_rate"]
