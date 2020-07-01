# IBM-Hackathon-2020 Summer Jam
Following is an attempt at a built-in package for surface codes to be pushed in the topological codes package of qiskit, by team Erwin's Tigers. The team memebers are @liuhenry @muirheadmaster @Phionx @shraggy and @zhenghaoding

Surface Code is a CSS code, consisting of pairwise commuting X and Z stabilizers made of Pauli gates. It creates a logical state on a 2 by 2 lattice made of quantum bits with the stabilizers X<sub>1</sub> X<sub>2</sub> Z<sub>1</sub> Z<sub>2</sub>. This repository has three parts: 
- **ciruits.py** creates initial circuit for measuring stabilizers and creating a logical state. Our code takes d (distance d) as input and T (no. of syndrome measurement rounds, usually T=d). 'd' should be odd and currently the code encodes only a logical 0 state. It's easy to make modifications and get logical 1,+,- states.
- **syndrome_graph** creates a graph which is a results of errors (as nodes) from varoius combination of pauli errors in the circuit for logical 0 from circuit.py. T
- **fitters.py** takes inputs from circuits.py and syndrome_graph.py and inserts errors in the original using a noise simulator. syndrome_graph.py is used to weigh the error graph to be sent in the minimum weight perfect matching algorithm to find disjoint edges and conclude a flip. The code returns the number of qubits flipped after all the stabilizer measurements, and concludes if there was a logical Z error in the final state.

## circuits.py
The SurfaceCode class, located in circuits.py creates the following circuit for any dstance, d. This example is for d=3 circuit where, blue patches are Z stabilizers and red patches are X stabilizers. Z stabilizers link qubits on data qubits (grey) in corners as control with syndrome qubits (black) in the centre.The Z and N marked in each patch determines the order of CX labeled in black. Grey coordinate labels are data qubit locations and black labels are syndrome qubit locations. The straight line marked Z<sub>L</sub> signifies that a logical Z is applied by operating Z on each qubit, on any horizontal line in the lattice. Similarly, the straight line marked X<sub>L</sub> signifies that a logical X is applied by operating X on each qubit, on any vertical line in the lattice. We choose one convention and say top edge signifies a Z logical operation and left edge signifies X logical operation!
<p align="center">
<img width="672" alt="Lattice" src="https://user-images.githubusercontent.com/293681/86267952-7541f700-bb95-11ea-8292-240bf344f7f8.png">
</p>
 The code creates the above circuit and measures each syndrome qubit. This is called syndrome measurement and is repeated T=d times. The results from each syndrome measurement are then processed to extract error nodes i.e. nodes which were flipped in consecutive syndrome measurements. This information is then utilised by the classes in syndrome_graph and fitter.py files to create error graphs and perform matching (explained in section "fitters.py"), to deduce the most probable error. Finally, logical Z error is concluded by checking if there were odd number of qubits with errors on top (Z<sub>L</sub>) edge and logical X error is concluded if there odd number of qubits with errors on the left (X<sub>L</sub>) edge
 
## Syndrome_Graph.py

The Syndrome class, located in syndrome_graph.py creates a graph of nodes from all the possible circuits by inserting an x error or a z error anywhere in the circuit. Following is an example of one such graph:
<p align="center">
<img width="361" alt="Graph" src="https://user-images.githubusercontent.com/293681/86267948-7410ca00-bb95-11ea-8c75-aacca29ceaa7.png">
</p>
Shortest path in this graph decides the weight of two edges when creating an "error graph" for GraphDecoder(in fitter.py). We analyse another method to do obtain a syndrom graph using graph traversal. This method is discussed in the next section.

## fitters.py
The GraphDecoder class, located in fitters.py, constructs the graph corresponding to the possible syndromes of a quantum error correction surface code, runs minimum weight perfect matching (MWPM) on a subgraph of detected errors to determine the most probable sequence of errors, and then flips a series of qubits to correct them. Our surface code is structured as a square lattice with alternating X and Z syndrome nodes, as depicted below for `d=5`:
<p align="center">
<img src="https://user-images.githubusercontent.com/42923017/86202361-01b9ce80-bb30-11ea-8656-820d8bb17085.jpg" width="50%">
</p>
This surface code evolves over a specified number of time steps `T`, effectively creating `T` syndrome node lattices. So, our 3D syndrome graph has the time step as the z dimension and the syndrome node lattices in the xy plane. Virtual nodes, which are nonphysical but nevertheless error-inducing, alternate between X and Z nodes across the border of the surface code and are also included in the graph. We construct our 3D syndrome graph by specifying the coordinates of the syndrome and virtual nodes, connecting the lattice of syndrome nodes at a given time step, and connecting the lattices between adjacent time steps and virtual nodes to their adjacent syndrome nodes at all time steps. Each edge weight is 1 to denote that adjacent syndrome nodes have a rotated Manhattan distance of 1. Below is an example of a 3D syndrome graph of X syndromes with `d=3` and `T=3`:
<p align="center">
<img src="https://user-images.githubusercontent.com/42923017/86195157-49375f00-bb1e-11ea-8dd1-63a3adae1002.jpg" width="50%">
</p>
Then, a subgraph of detected syndrome errors is extracted, where we account for path degeneracy in the edge weights and clone virtual nodes to allow for multiple virtual node to syndrome node matchings.
Below is an example of this error subgraph with `node_set = [(0, 1.5, 0.5), (1, 1.5, 0.5), (1, 0.5, 1.5), (2, 0.5, 1.5)]`:
<p align="center">
<img src="https://user-images.githubusercontent.com/42923017/86195162-4b99b900-bb1e-11ea-8f61-61ebf97a77f5.jpg" width="50%">
</p>
To determine the most probable set of syndrome errors, we run a MWPM on the error subgraph.
Below is an example of the MWPM matching graph for our error subgraph:
<p align="center">
<img src="https://user-images.githubusercontent.com/42923017/86195169-505e6d00-bb1e-11ea-9ad3-259d87911718.jpg" width="50%">
</p>
Finally, we correct the syndrome errors through a series of qubit flips.

# Acknowledgements
We would like to thank @quantumjim for valuable suggestions and pointing out to the topological_codes for repetition codes (also included here) repository availableon qiskit, which proved quite useful for our project.

# References
- [Surface Codes: Towards Practical Large-Scale Quantum Computation](https://arxiv.org/abs/1208.0928)
- [Stabilizer Codes and Quantum Error Correction](https://arxiv.org/pdf/quant-ph/9705052.pdf)
- [Multi-path Summation for Decoding 2D Topological Codes](https://quantum-journal.org/wp-content/uploads/2018/10/q-2018-10-19-102.pdf)
- [tutorial](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html#Lookup-table-decoding)
