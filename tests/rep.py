"""
Unittests for the rotated surface codes

python -m unittest tests/rep.py
"""
import sys
import unittest
from qiskit import execute, Aer
from qtcodes.circuits.base import LatticeError
from qtcodes.circuits.xxzz import XXZZQubit
from qtcodes.common import constants

from qtcodes.fitters.rotated_surface import RotatedDecoder
from tests.base import TestBase

sys.path.insert(0, "../")
from qtcodes import RepetitionQubit, RepetitionDecoder


class TestRep(TestBase, unittest.TestCase):
    """
    Test the Repetition Rotated Surface Code
    """

    encoder_type = XXZZQubit
    errors = ["x"]

    def setUp(self):
        self.params = {"d": 5}
        self.params["T"] = 1
        self.decoder = RepetitionDecoder(self.params)


class TestPhaseFlipProtectedRep(TestBase, unittest.TestCase):
    """
    Test the Repetition Rotated Surface Code
    """

    encoder_type = XXZZQubit
    errors = ["z"]

    def setUp(self):
        self.params = {"d": 5, "phase-flip-protected": True}
        self.params["T"] = 1
        self.decoder = RepetitionDecoder(self.params)


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

    def test_wrong_argument_type_phase_flip_protected(self):
        self.params = {}
        self.params["d"] = "1,3"
        self.params["phase-flip-protected"] = True
        self.params["T"] = 1

        with self.assertRaises(LatticeError):
            self.qubit = RepetitionQubit(self.params)

        with self.assertRaises(LatticeError):
            self.decoder = RepetitionDecoder(self.params)

    def test_wrong_height(self):
        self.params = {}
        self.params["d"] = (2, 3)
        self.params["phase-flip-protected"] = True
        self.params["T"] = 1

        with self.assertRaises(LatticeError):
            self.qubit = RepetitionQubit(self.params)

        with self.assertRaises(LatticeError):
            self.decoder = RepetitionDecoder(self.params)


# %%

if __name__ == "__main__":
    unittest.main()
