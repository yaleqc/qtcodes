# IBM-Hackathon-2020 summer jam
Following is an attempt at a built-in package for surface codes to be pushed in the topological codes package of qiskit.

The GraphDecoder class constructs the graph corresponding to the possible syndromes of a quantum error correction surface code, and then runs suitable decoders. The most probable errors are detected by a minimum weight perfect matching on the error graph, which consists of X or Z syndrome nodes and virtual nodes. When these errors are detected, they are corrected by a series of qubit flips.

# Sources

## Rep Code
- [tutorial](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html#Lookup-table-decoding)

