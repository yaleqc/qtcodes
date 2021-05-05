# Encoder Refactoring

## Directory structure

- `qiskit/ignis/verification/topological_codes`
  - `__init.py__`: import subdirectory classes and methods
  - `circuits`
    - `repetition_code.py`: formerly `circuits.py`
    - `surface_code.py`
  - `fitters`
    - `graph_decoder.py`: formerly `fitters.py`
    - `mwpm_decoder.py`
    - Misc methods in `fitters.py`

## `TopologicalCode` meta-class

### Abstract methods

#### General

- `__init__(d[, T])`
- `get_circuit_list()`
- `process_results(raw_results)`
- `readout()`

#### Circuit operations

- `syndrome_measurement([reset, basrrier])`
- `x([logs, barrier])`
- Other gates, need standardized name



