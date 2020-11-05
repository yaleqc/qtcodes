#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
simulation script for benchmark data
"""
#%%
import sys
import os

sys.path.insert(0, ".." + os.sep + ".." + os.sep)
from benchmarking.benchmarking_tools import SurfaceCodeBenchmarkingTool
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import pauli_error, depolarizing_error
from surface_code.fitters import GraphDecoder
from surface_code.circuits import SurfaceCodeLogicalQubit
from qiskit import QuantumCircuit, execute, QuantumRegister, ClassicalRegister, Aer
from tqdm import tqdm

import multiprocessing as mp

# Noise Model Function
def get_noise_model(p_err):

    error_gate1 = pauli_error([("X", p_err / 2), ("Z", p_err / 2), ("I", 1 - p_err)])

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_gate1, "id")
    return noise_model


if __name__ == "__main__":
    d = int(sys.argv[1])
    T = int(sys.argv[2])

    qubit = SurfaceCodeLogicalQubit(d)
    qubit.stabilize()
    qubit.identity_data()
    qubit.stabilize()
    qubit.readout_z()
    benchmarking_tool = SurfaceCodeBenchmarkingTool(
        decoder=GraphDecoder(d=d, T=T),
        readout_circuit=qubit,
        noise_model_func=get_noise_model,
    )

    print("\nSIMULATE: (d={},T={})\n".format(benchmarking_tool.d, benchmarking_tool.T))
    correct_logical_value = 0
    noise_values = [
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
    ]
    benchmarking_tool.simulate_readout_mp(
        correct_logical_value=correct_logical_value, noise_values=noise_values
    )

# %%
