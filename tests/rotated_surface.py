"""
Unit tests for the rotated surface codes

python -m unittest tests/rotated_surface.py
"""
from abc import ABCMeta
import unittest

from qtcodes.common import constants
from tests.base import TestBase


from qtcodes import XXZZQubit, XZZXQubit, RotatedDecoder


class TestSquareXXZZ(TestBase, unittest.TestCase):
    """
    Unit tests for the XXZZ (CSS) Rotated Surface Code
    """

    encoder_type = XXZZQubit

    def setUp(self):
        self.params = {"d": (5, 5)}
        self.params["T"] = 1
        self.decoder = RotatedDecoder(self.params)


class TestRectangularXXZZ(TestBase, unittest.TestCase):
    """
    Unit tests for the XXZZ (CSS) Rotated Surface Code
    """

    encoder_type = XXZZQubit

    def setUp(self):
        self.params = {"d": (3, 5)}
        self.params["T"] = 1
        self.decoder = RotatedDecoder(self.params)


class Test1DXXZZ(TestBase, unittest.TestCase):
    """
    Unit tests for the XXZZ (CSS) Rotated Surface Code
    """

    encoder_type = XXZZQubit

    def setUp(self):
        self.params = {"d": (5, 1)}
        self.params["T"] = 1
        self.decoder = RotatedDecoder(self.params)


class TestXZZX(TestBase, metaclass=ABCMeta):
    encoder_type = XZZXQubit

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
        dw = d[constants.DW]
        row = indx // dw
        col = indx % dw

        valid_syndrome = lambda x: self.decoder._valid_syndrome(
            x, "X"
        ) or self.decoder._valid_syndrome(x, "Z")

        if error_type == "x":
            neighbors = [
                (0.0, row - 0.5, col + 0.5),
                (0.0, row + 0.5, col - 0.5),
            ]
            return [x for x in neighbors if valid_syndrome(x[1:])]
        elif error_type == "z":
            neighbors = [
                (0.0, row - 0.5, col - 0.5),
                (0.0, row + 0.5, col + 0.5),
            ]
            return [x for x in neighbors if valid_syndrome(x[1:])]
        return []


class TestSquareXZZX(TestXZZX, unittest.TestCase):
    """
    Unit tests for the XZZX Rotated Surface Code
    """

    def setUp(self):
        self.params = {"d": (5, 5)}
        self.params["T"] = 1
        self.decoder = RotatedDecoder(self.params)


class TestRectangularXZZX(TestXZZX, unittest.TestCase):
    """
    Unit tests for the XZZX Rotated Surface Code
    """

    def setUp(self):
        self.params = {"d": (3, 5)}
        self.params["T"] = 1
        self.decoder = RotatedDecoder(self.params)


# %%

if __name__ == "__main__":
    unittest.main()
