# -*- coding: utf-8 -*-
"""
simulation script for benchmark data
"""
#%%
import sys
import os

sys.path.insert(0, ".." + os.sep)
from benchmarking.benchmarking_tools import SurfaceCodeBenchmarkingTool
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise.errors import pauli_error, depolarizing_error
from surface_code.fitters import GraphDecoder
from qiskit import QuantumCircuit, execute, QuantumRegister, ClassicalRegister, Aer
from tqdm import tqdm

import multiprocessing as mp

# Logical Readout Circuit (TODO REMOVE)
data = QuantumRegister(9, name="data")
mx = QuantumRegister(4, name="mx")
mz = QuantumRegister(4, name="mz")

rounds = 1  # The actual number of rounds will always be more than 1, since the first round creates the quiescent state
measurements = [
    ClassicalRegister(8, name="c{}".format(i + 1)) for i in range(rounds + 1)
]
logical = QuantumRegister(1, name="logical")

base_circ = QuantumCircuit(data, mz, mx, *measurements, logical)


def stabilize(circ, i):
    # Top left
    circ.h(mx[0])
    circ.cx(mx[0], data[1])
    circ.cx(mx[0], data[0])
    circ.cx(data[1], mz[0])
    circ.cx(data[0], mz[0])
    circ.cx(data[4], mz[0])
    circ.cx(data[3], mz[0])
    circ.h(mx[0])

    # Top right
    circ.h(mx[1])
    circ.cx(mx[1], data[2])
    circ.cx(mx[1], data[1])
    circ.cx(mx[1], data[5])
    circ.cx(mx[1], data[4])
    circ.cx(data[2], mz[1])
    circ.cx(data[5], mz[1])
    circ.h(mx[1])

    # Bottom left
    circ.h(mx[2])
    circ.cx(data[3], mz[2])
    circ.cx(data[6], mz[2])
    circ.cx(mx[2], data[4])
    circ.cx(mx[2], data[3])
    circ.cx(mx[2], data[7])
    circ.cx(mx[2], data[6])
    circ.h(mx[2])

    # Bottom right
    circ.h(mx[3])
    circ.cx(mx[3], data[8])
    circ.cx(mx[3], data[7])
    circ.cx(data[5], mz[3])
    circ.cx(data[4], mz[3])
    circ.cx(data[8], mz[3])
    circ.cx(data[7], mz[3])
    circ.h(mx[3])
    circ.barrier()

    circ.measure(mz, measurements[i][0:4])
    circ.measure(mx, measurements[i][4:8])
    circ.reset(mz)
    circ.reset(mx)
    circ.barrier()


def get_stabilized_circ(base_circuit, rounds):
    circ = base_circuit.copy()
    for i in range(rounds + 1):
        stabilize(circ, i)
    return circ


logical_zero = get_stabilized_circ(base_circ, rounds)
logical_readout = ClassicalRegister(1, name="logicalR")
logical_zero.add_register(logical_readout)

# The Z-logical readout
logical_zero.reset(logical)
logical_zero.cx(data[0], logical)
logical_zero.cx(data[1], logical)
logical_zero.cx(data[2], logical)
logical_zero.measure(logical, logical_readout)

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


decoder_keys = [(d, 1) for d in range(3, 13, 2)]
benchmarking_tools = []

for decoder_key in tqdm(decoder_keys):
    benchmarking_tools.append(
        SurfaceCodeBenchmarkingTool(
            decoder=GraphDecoder(d=decoder_key[0], T=decoder_key[1]),
            readout_circuit=logical_zero.copy(),
            noise_model_func=get_noise_model,
        )
    )


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
    simulate(i)

# Get mp working
# pool = mp.Pool(4)  # TODO change when running on HPC
# pool.map(test, list(range(4)))

# %%
