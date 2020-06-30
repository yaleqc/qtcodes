# IBM-Hackathon-2020 summer jam
Following is an attempt at a built-in package for surface codes to be pushed in the topological codes package of qiskit.

The GraphDecoder class, located in fitters.py, constructs the graph corresponding to the possible syndromes of a quantum error correction surface code, runs minimum weight perfect matching (MWPM) on a subgraph of detected errors to determine the most probable sequence of errors, and then flips a series of qubits to correct them. Our surface code is structured as a square lattice with alternating X and Z syndromes, as depicted below:
![Surface_Code_Diagram](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41534-018-0106-y/MediaObjects/41534_2018_106_Fig1_HTML.png?as=webp)

This surface code evolves over a specified number of time steps, giving rise to a 3D graph of X and Z syndrome nodes over time. Below is a diagram of how errors can occur over time across the surface code:
![Syndrome_3D_Graph](https://d3i71xaburhd42.cloudfront.net/c5cff89b63b34167235bd4ca2445b29c99ea07f8/4-Figure2-1.png)

We construct our 3D syndrome graph of either X or Z syndrome errors by specifying the coordinates of the syndrome nodes, connecting the lattice of syndrome nodes at a given time step, and connecting the lattices between adjacent time steps. Virtual nodes, which do not physically appear on the border of the diagram but can still lead to error, are also specified and connected to the appropriate border syndrome nodes. Then, a subgraph of detected errors is extracted, and we run a MWPM on the subgraph to determine the most probable errors. After these errors are detected, they are corrected by a series of qubit flips.

# Sources

## Rep Code
- [tutorial](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html#Lookup-table-decoding)

