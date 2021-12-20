from abc import ABCMeta
from qiskit import execute, Aer
from qtcodes.circuits.base import TQubit
from qtcodes.common import constants


class TestBase(metaclass=ABCMeta):
    """
    Abstract base class for testing
    """

    encoder_type: TQubit = None
    errors = ["x", "z"]

    def get_logical_error_rate(
        self,
        readout_strings,
        correct_logical_value,
        logical_readout_type,
        err_prob=None,
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
                readout, logical_readout_type, err_prob=err_prob
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
        dh = d[constants.DH]
        dw = d[constants.DW]
        row = indx // dw
        col = indx % dw
        if error_type == "x":
            parity = 2 * (indx % 2) - 1  # Even is -1, Odd is 1
            neighbors = [
                (0.0, row + parity * 0.5, col - 0.5),
                (0.0, row - parity * 0.5, col + 0.5),
            ]
            return [x for x in neighbors if x[1] > 0 and x[1] < dh - 1]
        elif error_type == "z":
            parity = 2 * (indx % 2) - 1  # Even is -1, Odd is 1
            neighbors = [
                (0.0, row - parity * 0.5, col - 0.5),
                (0.0, row + parity * 0.5, col + 0.5),
            ]
            return [x for x in neighbors if x[2] > 0 and x[2] < dw - 1]
        return []

    def simulate(self, error, i, logical_readout_type):
        # Set up circuit
        qubit = self.encoder_type(self.params)
        if logical_readout_type == "Z":
            qubit.reset_z()
        else:
            qubit.reset_x()
        qubit.stabilize()
        qubit.circ.__getattribute__(error)(qubit.lattice.qregisters["data"][i])  # error
        qubit.stabilize()
        if logical_readout_type == "Z":
            qubit.readout_z()
        else:
            qubit.readout_x()

        # Simulate
        results = (
            execute(qubit.circ, Aer.get_backend("aer_simulator"), shots=1000,)
            .result()
            .get_counts()
        )

        # Get Logical Error Rate
        logical_error_rate = self.get_logical_error_rate(
            results, 0, logical_readout_type
        )
        print(
            f"Decoding result for {error} Error on the {i}th data qubit ({logical_readout_type} readout), logical_error_rate: {logical_error_rate}."
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
        print(
            "\n"
            + 9 * "="
            + "\n"
            + self.__class__.__name__
            + ": "
            + str(self.params)
            + "\n"
        )
        # set up circuit
        d = self.params["d"]
        for i in range(d[constants.DH] * d[constants.DW]):
            for error in self.errors:
                if error == "x" and d[constants.DH] == 1:
                    continue

                if error == "z" and d[constants.DW] == 1:
                    continue

                self.simulate(error, i, "Z")
                self.simulate(error, i, "X")
