"""
TQubit constants
"""
from typing import Dict, Type, Any
from .xxzz import XXZZQubit
from .xzzx import XZZXQubit
from .repetition import RepetitionQubit

REPETITION = "Repetition"
XXZZ = "XXZZ"
XZZX = "XZZX"


blueprint: Dict[str, Type[Any]] = {
    REPETITION: RepetitionQubit,
    XXZZ: XXZZQubit,
    XZZX: XZZXQubit,
}
