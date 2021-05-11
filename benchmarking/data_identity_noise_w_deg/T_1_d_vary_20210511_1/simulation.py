# -*- coding: utf-8 -*-
"""
simulation script for benchmark data
"""
#%%
import sys
import os

sys.path.insert(0, ".." + os.sep + ".." + os.sep + ".." + os.sep)
from benchmarking.benchmarking_tools import SurfaceCodeBenchmarkingTool
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import pauli_error, depolarizing_error
from topological_codes import RotatedGraphDecoder
from topological_codes import XXZZQubit
from qiskit import QuantumCircuit, execute, QuantumRegister, ClassicalRegister, Aer
from tqdm import tqdm

import multiprocessing as mp

#%%
# Noise Model Function
def get_noise_model(p_err):

    error_gate1 = pauli_error([("X", p_err), ("I", 1 - p_err)])

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_gate1, "id")
    return noise_model


if __name__ == "__main__":
    decoder_keys = [(d, 1) for d in range(3, 11, 2)]
    benchmarking_tools = []

    for decoder_key in tqdm(decoder_keys):
        d = decoder_key[0]
        T = decoder_key[1]
        qubit = XXZZQubit({"d": d})
        qubit.stabilize()
        qubit.identity_data()
        qubit.stabilize()
        qubit.readout_z()
        benchmarking_tools.append(
            SurfaceCodeBenchmarkingTool(
                decoder=RotatedGraphDecoder({"d": d, "T": T}),
                readout_circuit=qubit.circ,
                noise_model_func=get_noise_model,
            )
        )
    print("\nDONE SETTING UP DECODERS!\n")

    for benchmarking_tool in tqdm(benchmarking_tools):
        print(
            "\nSIMULATE: (d={},T={})\n".format(benchmarking_tool.d, benchmarking_tool.T)
        )
        correct_logical_value = 0
        noise_values = [
            0.001,
            0.005,
            0.01,
            0.015,
            0.02,
            0.025,
            0.03,
        ]
        benchmarking_tool.simulate_readout_mp(
            correct_logical_value=correct_logical_value, noise_values=noise_values
        )

# %%
