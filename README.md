# qtcodes
*Qiskit Topological Codes*

[![License](https://img.shields.io/github/license/yaleqc/qtcodes.svg?style=popout-square)](https://opensource.org/licenses/Apache-2.0)
[![](https://img.shields.io/github/release/yaleqc/qtcodes.svg?style=popout-square)](https://github.com/yaleqc/qtcodes/releases)
[![](https://img.shields.io/pypi/dm/qtcodes.svg?style=popout-square)](https://pypi.org/project/qtcodes/)

## Installation

`qtcodes` is published on PyPI. So, to install, simply run:

```
pip install qtcodes
```

This will install a precompiled version of `qtcodes` into your python environment.



## Building from source

To build `qtcodes` from source, pip install using:

```
git clone https://github.com/yaleqc/qtcodes.git
cd qtcodes
pip install --upgrade .
```

To check if the installation was successful, run:
```
python3
>>> import qtcodes as qtc
```

---

## Motivation

Quantum computation is an inherently noisy process. Scalable quantum computers will require fault-tolerance to implement useful computation. There are many proposed approaches to this, but one promising candidate is the family of *topological quantum error correcting codes*.

Currently, the [`qiskit.ignis.verification.topological_codes`](https://qiskit.org/documentation/apidoc/verification.html#topological-codes) module provides a general framework for QEC and implements one specific example, the *repetition code*. Qiskit Topological Codes builds out the `qtcodes` module into a diverse family of QEC encoders and decoders, supporting the repetition code, XXXX/ZZZZ (XXZZ) rotated surface code, and the XZZX rotated surface code.

Inspired by the [Qiskit Textbook](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html), we've written a full set of [jupyter notebook tutorials](./tutorials) to demonstrate the [circuit encoders](./qtcodes/circuits), [graph decoders](./qtcodes/fitters), and [benchmarking tools](./qtcodes/tools/benchmarking.py) that compose Qiskit Topological Codes. These tutorials both demonstrate the elegance of QEC codes as well as the utility of this package -- please check them out!

## Codebase

<p align="center">
<img width="300" alt="surface code teaser" src="https://user-images.githubusercontent.com/10100490/130364540-58ec18b5-6e97-4625-a0ff-e8990f47782f.jpg"><br>
<div flush="left"><b>Fig 1.</b> Rotated XXXX/ZZZZ (XXZZ) Surface Code. ZZZZ/ZZ syndromes in red, XXXX/XX syndromes in purple, physical errors in green, and syndrome hits in yellow.</div>
</p>

Topological QEC codes disperse, and thus protect, one quantum bit of logical information across many physical qubits. The classical repetition code distributes 1 bit of logical information across multiple imperfect physical bits (e.g. logical 0 is 000...0 and logical 1 is 111...1). In the classical repetition logical 0 bit, for example, a few physical bits may flip to 1, but the majority will very likely stay in 0, thus preserving the logical 0 bit. Similarly, the surface code protects one logical qubit in a grid of imperfect physical qubits against Pauli errors.

The `qtcodes` module can be broken down into `circuits` (encoders) and `fitters` (decoders). Additionally, unittests can be found in `tests` and benchmarking tools in `qtcodes/tools`.


> The rotated surface code is based on the earlier theoretical idea of a [toric code](https://decodoku.blogspot.com/2016/03/6-toric-code.html), with periodic boundary conditions instead of open boundary conditions. This has been shown to be largely identical, but embedding a surface code on an actual device is much easier.

### Circuits

The `qtcodes.circuits` sub-module contains classes such as `XXZZQubit`, `XZZXQubit`, and `RepetitionQubit`, which each allow users to construct and manipulate circuits encoding one logical qubit using a particular QEC code.

For example, we can create and apply a logical X onto a `RepetitionQubit` as follows

```
from qtcodes import RepetitionQubit
qubit = RepetitionQubit({"d":3},"t")
qubit.reset_z()
qubit.stabilize()
qubit.x()
qubit.stabilize()
qubit.readout_z()
qubit.draw(output='mpl', fold=150)
```
![Repetition Code Qubit](https://user-images.githubusercontent.com/10100490/130364536-756c5e53-ef73-4564-8b44-6b978325fe9c.png)

`qtcodes.circuits.circ` also allows users to create  `TopologicalRegister`s (treg: a collection of topological qubits) and `TopologicalCircuit`s (tcirc: a circuit built using a treg), the analog of `QuantumRegister` and `QuantumCircuit`.

We can, for example, create a tcirc and treg out of two `RepetitionQubit`s.
```
from qtcodes import TopologicalRegister, TopologicalCircuit
treg = TopologicalRegister(2, ctype="Repetition", params={"d": 3})
circ = TopologicalCircuit(treg)
circ.x(treg[0])
circ.stabilize(treg[1])
circ.x(1)
circ.draw(output='mpl', fold=500)
```

![Repetition Code TCirc](https://user-images.githubusercontent.com/10100490/130364537-473de90a-6ed2-48a7-924f-00dd1fbce3e4.png)

Learn more about circuits through encoder tutorials such as this [one](./tutorials/xxzz/1-circuits.ipynb) for the XXXX/ZZZZ rotated surface code.

### Fitters

Topological codes aim to build better (read: less noisy) logical qubits out of many imperfect physical qubits. This improvement is enabled by decoding schemes that can detect and thus correct for errors on a code's constituent physical qubits.

The Qiskit Topological Codes package leverages Minimum-Weight Perfect Matching Graph Decoding to efficiently correct logical qubit readout.


For example, we can decode the syndrome hits in Fig 1 and fine the most probable error chains (data qubit flips) corresponding to these syndrome hits.
```
#d: surface code side length, T: number of rounds
decoder = RotatedDecoder({"d":5,"T":1})
all_syndromes = {"X": [(0,1.5,.5),(0,.5,1.5)], "Z": [(0,0.5,0.5),(0,1.5,1.5),(0,1.5,3.5), (0,3.5,3.5)]}
matches = {}

for syndrome_key, syndromes in all_syndromes.items():
    print(f"{syndrome_key} Syndrome Graph")
    error_graph = decoder._make_error_graph(syndromes,syndrome_key)
    print("Error Graph")
    decoder.draw(error_graph)
    matches[syndrome_key] = decoder._run_mwpm(error_graph)
    matched_graph = decoder._run_mwpm_graph(error_graph)
    print("Matched Graph")
    decoder.draw(matched_graph)
    print(f"Matches: {matches[syndrome_key]}")
    print("\n===\n")
```

<p align="middle" style="background:#fff">
  <img src="https://user-images.githubusercontent.com/10100490/130364532-93f60a0f-1636-4967-b324-6745e23a003a.png" width="49%" align="top"/>
  <img src="https://user-images.githubusercontent.com/10100490/130364534-4316ef5a-38f4-4d67-83b9-dde62e2bf8c2.png" width="49%" align="top"/>
</p>

In this way, Qiskit Topological Codes uses graph decoding to find and correct for the most probable set of errors (error chains).

The careful reader will notice that connecting syndrome hits in the most probable set of "error chains" does not uniquely specify the underlying physical qubits that underwent physical errors (i.e. there are multiple shortest paths between two syndrome hits). It turns out, by the nuances of how topological codes store  logical information (i.e. codespace), in most cases the exact path across physical qubits doesn't matter when correcting for an error chain. Read more about this in this [tutorial](./tutorials/xxzz/2-fitters.ipynb) on Decoding for XXZZ Qubits!

### Benchmarking

Finally, the efficiency and efficacy of the Qiskit Topological Codes package is demonstrated through benchmark simulations achieving threshold for the Repetition, XXZZ, and XZZX topological codes. Here, threshold is defined as the maximum physical error rate (i.e. imperfection level of physical qubits) below which larger surface codes perform better than smaller surface codes.

<p align="middle">
  <img src="https://user-images.githubusercontent.com/10100490/130364554-7e1536f2-a6be-4487-b0dc-651af1a861b9.png" width="32%" />
  <img src="https://user-images.githubusercontent.com/10100490/130364555-3f84d008-851b-42db-a8fc-e9af5864586d.png" width="32%" />
  <img src="https://user-images.githubusercontent.com/10100490/130364556-32375c1c-2be0-4707-a44d-cfa080265800.png" width="32%" /><br>
  <div flush="left">
  <b>Fig. 2</b> By simulating circuits with errors inserted between two rounds of stabilizing measurements, we are able to extract a logical error rate for each code for a given physical error rate (quality of physical qubit) and surface code size. In particular, threshold is shown for the repetition code (left), XXZZ code (center), and XZZX code (right).</div>
</p>

Explore the benchmarking [tools](./qtcodes/tools/benchmarking.py) and [simulations](./data/) to see how the graphs in Fig. 2 were created.

## Future Directions

*Checkout [issues](https://github.com/yaleqc/qtcodes/issues) to see what we are working on these days!*

## Acknowledgements

**Core Devs:** [Shantanu Jha](https://github.com/Phionx), [Jessie Chen](https://github.com/JazzyCH), [Aaron Householder](https://github.com/aaronhouseholder), [Allen Mi](https://github.com/Allenator)

Thanks to our mentor [James Wootton](https://github.com/quantumjim) (IBM) for invaluable feedback and support since the inception of this project at the IBM Qiskit - Summer Jam Hackathon 2020.

Thanks also to [Matthew Treinish](https://github.com/mtreinish) from the [retworkx](https://github.com/Qiskit/retworkx) team for helping onboard and support this project.

**Alums:** [Henry Liu](https://github.com/liuhenry), [Shraddha Singh](https://github.com/shraggy), [Will Sun](https://github.com/muirheadmaster), [Andy Ding](https://github.com/ZhenghaoDing)


## References

*Here's some reading material that we found particularly useful:*
- Presentation [slides](https://docs.google.com/presentation/d/1HC5tQkvOcfl5lPDWy-l8ZvsWW3AGf6jUCmMu_sz-cDY/edit?usp=sharing) and [video](https://www.youtube.com/watch?v=jb1qD0pZbF4&list=PLOFEBzvs-VvqQMAVaXoFlSqjqgbX5k-fL&index=18&ab_channel=QiskitQiskit) about this package to the Qiskit Advocate community at the November 2020 Qiskit Advocate Demo Session.
- [Surface Codes: Towards Practical Large-Scale Quantum Computation](https://arxiv.org/abs/1208.0928)
- [Stabilizer Codes and Quantum Error Correction](https://arxiv.org/pdf/quant-ph/9705052.pdf)
- [Multi-path Summation for Decoding 2D Topological Codes](https://quantum-journal.org/wp-content/uploads/2018/10/q-2018-10-19-102.pdf)
- [Qiskit Textbook - Introduction to Quantum Error Correction using Repetition Codes](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html)
