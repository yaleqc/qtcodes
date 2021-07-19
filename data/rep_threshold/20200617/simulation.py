"""
Simulation script to demonstrate Rep Code Threshold
"""

import sys
import os
from tqdm import tqdm
from qiskit.providers.aer.noise.errors import pauli_error
from qiskit.providers.aer.noise import NoiseModel

sys.path.insert(0, ".." + os.sep + ".." + os.sep + ".." + os.sep)
from qtcodes import TopologicalBenchmark
from qtcodes import RepetitionDecoder, RepetitionQubit


# Noise Model Function
def noise_model_func(p_err):

    error_gate1 = pauli_error([("X", p_err), ("I", 1 - p_err)])

    noise_model = NoiseModel()
    noise_model.add_all_qubit_quantum_error(error_gate1, "id")
    return noise_model


if __name__ == "__main__":
    decoder_keys = [(d, 1) for d in range(3, 13, 2)]
    tools = []
    print("Setting up decoders!\n")
    for decoder_key in tqdm(decoder_keys):
        d = decoder_key[0]
        T = decoder_key[1]

        # setup circ
        qubit = RepetitionQubit({"d": d})
        qubit.reset_z()
        qubit.stabilize()
        qubit.id_data()
        qubit.stabilize()
        qubit.readout_z()

        tools.append(
            TopologicalBenchmark(
                decoder=RepetitionDecoder({"d": d, "T": T}),
                circ=qubit.circ,
                noise_model_func=noise_model_func,
                correct_logical_value=0,
            )
        )

    for tool in tools:
        print(
            "\nSIMULATE: (d={},T={})\n".format(
                tool.decoder.params["d"], tool.decoder.params["T"]
            )
        )
        physical_error_rates = [
            0.05,
            0.15,
            0.25,
            0.35,
            0.45,
            0.55,
            0.65,
            0.75,
            0.85,
            0.95,
        ]
        tool.sweep(
            physical_error_rates=physical_error_rates, shots=1024 * 16,
        )
