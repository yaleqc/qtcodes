# -*- coding: utf-8 -*-
"""
benchmarking class for qtcodes
"""
import sys
import os
import glob

sys.path.insert(0, ".." + os.sep)

import numpy as np
import matplotlib.pyplot as plt

from qiskit import execute, Aer, QuantumCircuit
from tqdm import tqdm

from typing import Dict, Optional, List


plt.rcParams.update({"font.size": 14, "pdf.fonttype": 42, "ps.fonttype": 42})

from multiprocessing import Pool


class TopologicalBenchmark:
    def __init__(
        self,
        decoder,
        circ: QuantumCircuit,
        noise_model_func=None,
        correct_logical_value: int = 0,
    ):
        self.decoder = decoder
        self.circ = circ
        self.filename = f"d_{self.decoder.params['d']}_T_{self.decoder.params['T']}.npz"

        self.correct_logical_value = correct_logical_value
        self.noise_model_func = noise_model_func
        self.data: Dict[str, List[float]] = {
            "physical_error_rates": [],
            "logical_error_rates": [],
        }

    def logical_error_rate(
        self, readout_strings: Dict[str, int], err_prob: Optional[float] = None
    ) -> float:
        """
        Args:
            readout_strings: a dictionary of readout strings along with counts
            e.g. {"1 00000000 00000000":48, "1 00100000 00100000":12, ...} in the case of d=3 and T=2

            err_prob (Optional[float]): Probability of IID data qubit X/Z flip. Defaults to None.

        Returns:
            error_rate (float): = (number of unsuccessful logical value predictions) / (total number of predictions )
        """
        total_count = 0.0
        total_errors = 0.0
        for readout, count in readout_strings.items():
            total_count += count
            predicted_logical_value = self.decoder.correct_readout(
                readout, "Z", err_prob=err_prob
            )
            if predicted_logical_value != self.correct_logical_value:
                total_errors += count

        return total_errors / total_count

    def sweep(
        self,
        physical_error_rates: Optional[List[float]] = None,
        save_data: bool = True,
        shots: int = 2048,
        deg_weight: bool = True,
    ) -> None:
        """
        Sweep physical error rates and calculate the associated logical error rate.

        Args:
            physical_error_rates (Optional[List[float]]):
                List of physical error rates to sweep.

            save_data (bool):
                This boolean determines whether this data is saved to an npz file.

            shots (int):
                Shots in the circuit simulation.

            deg_weight (bool):
                Whether or not to use degeneracy weighting.

        """
        self.data["physical_error_rates"] = []
        self.data["logical_error_rates"] = []

        physical_error_rates = sorted(
            physical_error_rates
            if physical_error_rates is not None
            else [0.04, 0.07, 0.10, 0.13, 0.16,],
            reverse=True,
        )  # higher physical_error_rate readout is slower to decode, gives more accurate tqdm estimate
        pbar = tqdm(physical_error_rates)
        for physical_error_rate in pbar:
            results = (
                execute(
                    self.circ,
                    Aer.get_backend("aer_simulator"),
                    noise_model=self.noise_model_func(physical_error_rate),
                    shots=shots,
                )
                .result()
                .get_counts()
            )
            logical_error_rate_value = self.logical_error_rate(
                results, err_prob=physical_error_rate if deg_weight else None
            )
            self.data["physical_error_rates"].append(physical_error_rate)
            self.data["logical_error_rates"].append(logical_error_rate_value)
            if save_data:
                self.append_data(physical_error_rate, logical_error_rate_value)
            pbar.set_description(f"Done with noise: {physical_error_rate}")

        # final sort
        physical_error_rates_final = np.array(self.data["physical_error_rates"])
        logical_error_rates_final = np.array(self.data["logical_error_rates"])
        indxs = np.argsort(physical_error_rates_final)
        self.data["physical_error_rates"] = physical_error_rates_final[indxs]
        self.data["logical_error_rates"] = logical_error_rates_final[indxs]

    def sweep_mp(
        self,
        physical_error_rates: Optional[List[float]] = None,
        save_data: bool = True,
        shots: int = 2048,
    ) -> None:
        """
        Multi-processed weep physical error rates and calculate the associated logical error rate.

        Args:
            physical_error_rates (Optional[List[float]]):
                List of physical error rates to sweep.

            save_data (bool):
                This boolean determines whether this data is saved to an npz file.

            shots (int):
                Shots in the circuit simulation.

        """
        physical_error_rates = sorted(
            physical_error_rates
            if physical_error_rates is not None
            else [0.04, 0.07, 0.10, 0.13, 0.16,],
        )
        p = Pool(4)  # TODO change on HPC
        single_lambda = lambda physical_error_rate: self.single(
            physical_error_rate, save_data=save_data, shots=shots
        )
        logical_error_rates = p.map(single_lambda, physical_error_rates)

        self.data["physical_error_rates"] = physical_error_rates
        self.data["logical_error_rates"] = logical_error_rates

    def single(
        self, physical_error_rate: float, save_data: bool = True, shots: int = 2048,
    ):
        """
        Take single error rates and calculate the associated logical error rate.

        Args:
            physical_error_rate (float):
                Single physical error rate

            save_data (bool):
                This boolean determines whether this data is saved to an npz file.

            shots (int):
                Shots in the circuit simulation.

        """
        results = (
            execute(
                self.circ,
                Aer.get_backend("aer_simulator"),
                noise_model=self.noise_model_func(physical_error_rate),
                shots=shots,
            )
            .result()
            .get_counts()
        )
        logical_error_rate_value = self.logical_error_rate(
            results, err_prob=physical_error_rate
        )
        print("Done simulating physical_error_rate: " + str(physical_error_rate))
        if save_data:
            self.append_data(physical_error_rate, logical_error_rate_value)
        return logical_error_rate_value

    def append_data(self, physical_error_rate: float, logical_error_rate: float):
        try:
            data = np.load(self.filename)
            physical_error_rates = data["physical_error_rates"]
            logical_error_rates = data["logical_error_rates"]
        except:
            physical_error_rates = np.array([])
            logical_error_rates = np.array([])

        physical_error_rates = np.append(physical_error_rates, physical_error_rate)
        logical_error_rates = np.append(logical_error_rates, logical_error_rate)

        indxs = np.argsort(physical_error_rates)
        np.savez(
            self.filename,
            d=self.decoder.params["d"],
            T=self.decoder.params["T"],
            physical_error_rates=physical_error_rates[indxs],
            logical_error_rates=logical_error_rates[indxs],
        )


class TopologicalAnalysis:
    def __init__(self, filename: Optional[str] = None):
        self.filename = filename
        self.data: Dict[str, List[float]] = {}
        self.params: Dict[str, int] = {}

    def load_data(self):
        data = np.load(self.filename)
        self.params["d"] = int(data["d"])
        self.params["T"] = int(data["T"])
        self.data["physical_error_rates"] = data["physical_error_rates"]  # TODO
        self.data["logical_error_rates"] = data["logical_error_rates"]  # TODO

    def plot(self, fig=None, ax=None, log=True, per_round=False, **kwargs):
        if fig is None:
            fig = plt.figure(figsize=(3.5, 2.5), dpi=200)
        if ax is None:
            ax = fig.subplots()

        plt.plot(
            np.array(self.data["physical_error_rates"]),
            np.array(self.data["logical_error_rates"])
            / (self.params["T"] if per_round else 1.0),
            ".:",
            **kwargs,
        )
        if log:
            plt.yscale("log")
            plt.xscale("log")


class TopologicalBatchAnalysis:
    def __init__(self, dirname: str):
        self.dirname = dirname
        self.filenames = glob.glob(self.dirname + "*.npz")
        self.analyses = []
        for filename in self.filenames:
            self.analyses.append(TopologicalAnalysis(filename))
            self.analyses[-1].load_data()

        # sort analysis by increasing "d"
        sorted_indxs = np.argsort(
            np.array([analysis.params["d"] for analysis in self.analyses])
        )
        sorted_analyses = [self.analyses[i] for i in sorted_indxs]
        self.analyses = sorted_analyses

    def plot(self):
        for log_plot in [True, False]:
            # Plotting
            fig = plt.figure(figsize=(3.5, 2.5), dpi=200)
            ax = fig.subplots()

            for analysis in self.analyses:
                analysis.plot(
                    fig=fig,
                    ax=ax,
                    log=log_plot,
                    label="d={},T={}".format(
                        analysis.params["d"], analysis.params["T"]
                    ),
                )

            plt.plot(
                self.analyses[0].data["physical_error_rates"],
                self.analyses[0].data["physical_error_rates"],
                "-",
                label="breakeven",
            )
            plt.legend(loc="upper left", prop={"size": 6})
            ax.set_xlabel("Physical Error Rate", size=10)
            ax.set_ylabel("Logical Error Rate", size=10)
            ax.set_title("Comparison of Codes", size=10)
            fig.tight_layout()
            plt.savefig(
                self.dirname + "comparison" + ("_log" if log_plot else "") + ".png"
            )
            plt.show()

