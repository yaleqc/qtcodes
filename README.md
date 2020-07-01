# IBM-Hackathon-2020 Summer Jam
Following is an attempt at a built-in package for surface codes to be pushed in the topological codes package of qiskit.

The GraphDecoder class, located in fitters.py, constructs the graph corresponding to the possible syndromes of a quantum error correction surface code, runs minimum weight perfect matching (MWPM) on a subgraph of detected errors to determine the most probable sequence of errors, and then flips a series of qubits to correct them. Our surface code is structured as a square lattice with alternating X and Z syndrome nodes, as depicted below for `d=3` and `d=5`:

<img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41534-018-0106-y/MediaObjects/41534_2018_106_Fig1_HTML.png?as=webp" width="50%">

This surface code evolves over a specified number of time steps `T`, giving rise to a 3D syndrome graph of the syndrome node lattice over time. Virtual nodes, which are nonphysical but nevertheless error-inducing, alternate between X and Z nodes across the border of the surface code and are also included in the graph. We construct our 3D syndrome graph by specifying the coordinates of the syndrome and virtual nodes, connecting the lattice of syndrome nodes at a given time step, and connecting the lattices between adjacent time steps and virtual nodes to their adjacent syndrome nodes at all time steps. Each edge weight is 1 to denote that adjacent syndrome nodes have a rotated Manhattan distance of 1. Below is an example of a 3D syndrome graph of X syndromes with `d=3` and `T=3`:

<img src="https://user-images.githubusercontent.com/42923017/86195157-49375f00-bb1e-11ea-8dd1-63a3adae1002.jpg" width="50%">

Then, a subgraph of detected syndrome errors is extracted, where we account for path degeneracy in the edge weights and clone virtual nodes to allow for multiple virtual node to syndrome node matchings.
Below is an example of this error subgraph with `node_set = [(0, 1.5, 0.5), (1, 1.5, 0.5), (1, 0.5, 1.5), (2, 0.5, 1.5)]`:

<img src="https://user-images.githubusercontent.com/42923017/86195162-4b99b900-bb1e-11ea-8f61-61ebf97a77f5.jpg" width="50%">

To determine the most probable set of syndrome errors, we run a MWPM on the error subgraph.
Below is an example of the MWPM matching graph for our error subgraph:

<img src="https://user-images.githubusercontent.com/42923017/86195169-505e6d00-bb1e-11ea-9ad3-259d87911718.jpg" width="50%">

Finally, we correct the syndrome errors through a series of qubit flips.

# Sources
- [Surface Codes: Towards Practical Large-Scale Quantum Computation](https://arxiv.org/abs/1208.0928)
- [Stabilizer Codes and Quantum Error Correction](https://arxiv.org/pdf/quant-ph/9705052.pdf)
- [Multi-path Summation for Decoding 2D Topological Codes](https://quantum-journal.org/wp-content/uploads/2018/10/q-2018-10-19-102.pdf)
## Rep Code
- [tutorial](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html#Lookup-table-decoding)

