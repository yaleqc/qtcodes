"""
Unittests for the rotated surface codes
"""
import sys
import unittest
from qiskit import execute, Aer

sys.path.insert(0, "../")
from topological_codes import XXZZQubit, RotatedDecoder


class TestXXZZ(unittest.TestCase):
    """
    Unittests for the XXZZ (CSS) Rotated Surface Code
    """

    def setUp(self):
        self.params = {"d": 5}
        self.params["T"] = 1
        self.decoder = RotatedDecoder(self.params)

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

    def get_neighbors(self, indx: int, error_type: str):
        """
        Returns the syndrome node positions given some data qubit index
        and error_type on that data qubit.

        Args:
            indx (int): index of data qubit
            error_type (str): either "x" or "z" error on data qubit

        Returns:
            neighbors (List[Tuple[float]]):
                List of neighboring syndrome nodes
                that would be set off by the specified error
                on the specified data qubit.
        """
        d = self.params["d"]
        row = indx // d
        col = indx % d
        if error_type == "x":
            parity = 2 * (indx % 2) - 1  # Even is -1, Odd is 1
            neighbors = [
                (0.0, row + parity * 0.5, col - 0.5),
                (0.0, row - parity * 0.5, col + 0.5),
            ]
            return [x for x in neighbors if x[1] > 0 and x[1] < d - 1]
        elif error_type == "z":
            parity = 2 * (indx % 2) - 1  # Even is -1, Odd is 1
            neighbors = [
                (0.0, row - parity * 0.5, col - 0.5),
                (0.0, row + parity * 0.5, col + 0.5),
            ]
            return [x for x in neighbors if x[2] > 0 and x[2] < d - 1]
        return []

    def test_single_errors(self):
        """
        Setting up a |+z> state, stabilizing twice, and reading out logical Z.
        Inserting X and Z gates on each data qubit between tbe two stabilizer measurement rounds.

        Then, testing:
        1.  The MWPM decoder is able to correct these single qubit errors
            and predict the correct logical readout value.
        2.  Checking if the _string2node readout string to
            syndrome node parser is working correctly.
        """
        # set up circuit
        for i in range(self.params["d"] ** 2):
            for error in ["x", "z"]:
                # Set up circuit
                qubit = XXZZQubit(self.params)
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

                # Check if syndromes were parsed correctly
                expected_neighbors = self.get_neighbors(i, error)
                print(f"Expected neighbors: {expected_neighbors}")

                readout_strings = list(results.keys())
                for readout_string in readout_strings:
                    _, processed_results = self.decoder.parse_readout(readout_string)
                    nodes = sum(list(processed_results.values()), [])
                    for x in nodes:
                        self.assertIn(x, expected_neighbors)


# %%

if __name__ == "__main__":
    unittest.main()
