"""
Unittests for the rotated surface codes

python -m unittest tests/rep.py
"""
import sys
import unittest
from qiskit import execute, Aer
from qtcodes.circuits.base import LatticeError
from qtcodes.common import constants

from qtcodes.fitters.rotated_surface import RotatedDecoder

sys.path.insert(0, "../")
from qtcodes import RepetitionQubit, RepetitionDecoder


class TestRep(unittest.TestCase):
    """
    Test the Repetition Rotated Surface Code
    """

    def setUp(self):
        self.params = {"d": 5}
        self.params["T"] = 1
        self.decoder = RepetitionDecoder(self.params)

    def get_logical_error_rate(
        self, readout_strings, correct_logical_value, err_prob=None
    ):
        """
        Args:
            readout_strings: a dictionary of readout strings along with counts
            e.g. {"1 00000000 00000000":48, "1 00100000 00100000":12, ...} in the case of d=3, T=2

            correct_logical_value: integer (0/1) depicting original encoded logical value

        Returns:
            error_rate:
                float = (# of unsuccessful logical value predictions) / (total # of pred )
        """
        total_count = 0
        total_errors = 0
        for readout, count in readout_strings.items():
            total_count += count
            predicted_logical_value = self.decoder.correct_readout(
                readout, "Z", err_prob=err_prob
            )
            if predicted_logical_value != correct_logical_value:
                total_errors += count

        return total_errors / total_count

    def test_single_errors_rep(self):
        """
        Setting up a |+z> state, stabilizing twice, and reading out logical Z.
        Inserting X gates on each data qubit between tbe two stabilizer measurement rounds.

        Then, testing:
        1.  The MWPM decoder is able to correct these single qubit errors
            and predict the correct logical readout value.
        """
        # set up circuit
        d = self.params["d"]
        for i in range(d[constants.DH] * d[constants.DW]):
            for error in ["x"]:
                # Set up circuit
                qubit = RepetitionQubit(self.params)
                qubit.reset_z()
                qubit.stabilize()
                qubit.circ.__getattribute__(error)(
                    qubit.lattice.qregisters["data"][i]
                )  # error
                qubit.stabilize()
                qubit.readout_z()

                # Simulate
                results = (
                    execute(qubit.circ, Aer.get_backend("aer_simulator"), shots=1000,)
                    .result()
                    .get_counts()
                )
                print(results)

                # Get Logical Error Rate
                logical_error_rate = self.get_logical_error_rate(results, 0)
                print(
                    f"Decoding result for {error} Error on the {i}th data qubit, logical_error_rate: {logical_error_rate}."
                )
                self.assertEqual(
                    logical_error_rate,
                    0,
                    f"Decoding did not work for an {error} Error on the {i}th data qubit.",
                )


class TestRepExceptions(unittest.TestCase):
    """
    Test the Repetition Rotated Surface Code exceptions if wrong arguments are passed
    """

    def test_wrong_argument_type(self):
        self.params = {}
        self.params["d"] = "3,1"
        self.params["T"] = 1

        with self.assertRaises(LatticeError):
            self.qubit = RepetitionQubit(self.params)

        with self.assertRaises(LatticeError):
            self.decoder = RepetitionDecoder(self.params)

    def test_wrong_width(self):
        self.params = {}
        self.params["d"] = (3, 2)
        self.params["T"] = 1

        with self.assertRaises(LatticeError):
            self.qubit = RepetitionQubit(self.params)

        with self.assertRaises(LatticeError):
            self.decoder = RepetitionDecoder(self.params)


# %%

if __name__ == "__main__":
    unittest.main()
