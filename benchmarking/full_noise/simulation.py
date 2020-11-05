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
    error_gate2 = error_gate1.tensor(error_gate1)

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(
        error_gate1, "measure"
    )  # measurement error is applied to measurements
    noise_model.add_all_qubit_quantum_error(
        error_gate1, ["u1", "u2", "u3"]
    )  # single qubit gate error is applied to x gates
    noise_model.add_all_qubit_quantum_error(
        error_gate2, ["cx"]
    )  # two qubit gate error is applied to cx gates

    return noise_model


decoder_keys = [(d, 1) for d in range(3, 11, 2)]
benchmarking_tools = []

for decoder_key in tqdm(decoder_keys):
    d = decoder_key[0]
    T = decoder_key[1]
    qubit = SurfaceCodeLogicalQubit(d)
    qubit.stabilize()
    qubit.stabilize()
    qubit.readout_z()
    benchmarking_tools.append(
        SurfaceCodeBenchmarkingTool(
            decoder=GraphDecoder(d=d, T=T),
            readout_circuit=qubit,
            noise_model_func=get_noise_model,
        )
    )
print("\nDONE SETTING UP DECODERS!\n")


def simulate(index):
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
    benchmarking_tools[index].simulate_readout(
        correct_logical_value=correct_logical_value, noise_values=noise_values
    )


for i in range(len(benchmarking_tools)):
    print(
        "\nSIMULATE: (d={},T={})\n".format(
            benchmarking_tools[i].d, benchmarking_tools[i].T
        )
    )
    simulate(i)

# Get mp working
# pool = mp.Pool(4)  # TODO change when running on HPC
# pool.map(test, list(range(4)))

# %%
