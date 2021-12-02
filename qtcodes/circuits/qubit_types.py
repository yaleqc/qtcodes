"""
Qubit Types
"""
from typing import Dict, Type, Any
from qtcodes.circuits.xxzz import XXZZQubit
from qtcodes.circuits.xzzx import XZZXQubit
from qtcodes.circuits.repetition import RepetitionQubit
from qtcodes.common.constants import REPETITION, XZZX, XXZZ

str2qtype: Dict[str, Type[Any]] = {
    REPETITION: RepetitionQubit,
    XXZZ: XXZZQubit,
    XZZX: XZZXQubit,
}

