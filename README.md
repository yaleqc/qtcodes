# Qiskit Surface Code Encoder/Decoder
## IBM Qiskit - Summer Jam Hackathon 2020

Quantum computation is an inherently noisy process. Scalable quantum computers will require fault-tolerance to implement useful computation. There are many proposed approaches to this, but one promising candidate is the family of *topological quantum error correcting codes*.

Currently, the [`qiskit.ignis.verification.topological_codes`](https://qiskit.org/documentation/apidoc/verification.html#topological-codes) module provides a general framework for QEC and implements one specific example, the *repetition code*.

For the hackathon, our team [Erwin's Tigers](#team) implemented a **surface code encoder and decoder** for Qiskit Ignis. We hope that this implementation will be useful to other Qiskitters and will inspire others to continue building out the `topological_codes` module into a diverse family.

## Usage
Inspired by the [Qiskit Textbook](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html), we've written a full set of [jupyter notebook tutorials](https://github.com/Phionx/qiskit_surface_codes/tree/master/tutorials), which are the best way to get up to speed. They detail both the API and the gritty implementation details -- please check them out!

## Background

<p align="center">
<img width="300" alt="surface code teaser" src="https://user-images.githubusercontent.com/293681/86277823-5f3c3280-bba5-11ea-9c87-d8525a8cbe88.jpg">
</p>

Surface codes are a type of CSS code, consisting of pairwise commuting X and Z stabilizers made of Pauli gates. It defines a logical state on a 2 by 2 lattice made of quantum bits with the stabilizers X<sub>1</sub> X<sub>2</sub> Z<sub>1</sub> Z<sub>2</sub>.

The code is based on the earlier theoretical idea of a *toric code*, with periodic boundary conditions instead of open boundary conditions. This has been shown to be largely identical, but embedding a surface code on an actual device is much easier.

## Implementation

In general, we try to follow the existing structure of [`qiskit.ignis.verification.topological_codes`](https://qiskit.org/documentation/apidoc/verification.html#topological-codes). The code is implemented separately here but is able to easily be merged into Ignis.

There are two main interfaces — corresponding to the encoder and decoder, respectively:

### `SurfaceCode` in [`surface_code.circuits`](surface_code/circuits.py)

`SurfaceCode(d, T)` generates a `QuantumCircuit` for creating a logical state and measuring stabilizers. The class is parameterized with the code distance `d` (which should be odd) and the number of syndrome measurement rounds `T` (usually `T = d`). This class also handles parsing of the physical device readout into a form suitable for decoding. Please see the [encoder tutorial](tutorials/1_surface_code_encoding.ipynb) for a full walkthrough.
<p align="center">
<img width="615" alt="circuit" src="https://user-images.githubusercontent.com/293681/86277098-23ed3400-bba4-11ea-8305-c6d19eb73899.png">
</p>

### `GraphDecoder` in [`surface_code.fitters`](surface_code/fitters.py)

`GraphDecoder(d, T)` implements minimum-weight perfect matching (MWPM) on the syndrome measurements of the physical circuit. This class is similar to the existing `GraphDecoder` for repetition codes, but introduces a new framework to handle the 2D lattice.

Parsed readout from the device is used to generate graphs of *error chains* in time and space, which decode syndrome measurements into the most likely sequence of qubit flips over time. Please see the [decoder tutorial](tutorials/2_surface_code_decoding.ipynb) for a full walkthrough.

<p align="center">
<img width="615" alt="matching graph" src="https://user-images.githubusercontent.com/293681/86277350-8ba37f00-bba4-11ea-9560-02d5ea3167cd.png">
</p>

## Future Directions
The scope of the project is quite large, so we focused on completing a "minimum viable product" during the hackathon. However, there are many areas which we'd like to explore going forward. A few immediate ones:

* Expand `SurfaceCode(d, T).circuits` a full set of logical states (1, +, -) -- and ultimately logical gate operations for computation.
* Full benchmark of the physical-error to logical-error probabilities to determine the error correction threshold
* More simulation runs: different X/Z error probabilities, more limited `coupling_map`, etc.
* Our MWPM matching already has the below improvements to the basic algorithm, but are there more?
  * For a given pair of syndromes, there may be many possible error chains through space and time. We compute this "path degeneracy" and use it to re-weight the error probabilities.
  * We cross-match X and Z errors to produce an overall Y error. However, this doesn't exactly match a depolarizing channel, so ideally the weights would be re-adjusted with conditional probabilities.
* Other approaches to error-chain matching (e.g. [neural networks](https://iopscience.iop.org/article/10.1088/2058-9565/aa955a/meta) or [tensor networks](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.90.032326))?
* Our `GraphDecoder` implements two different approaches to syndrome graph generation. One is a "analytic" approach (much faster), and the other uses simulation to insert errors into the circuit. These produce slightly different syndrome graphs, but we get the same decoding results in our tests.

## Team
* [Andy Ding](https://github.com/ZhenghaoDing)
* [Shantanu Jha](https://github.com/Phionx)
* [Henry Liu](https://github.com/liuhenry)
* [Shraddha Singh](https://github.com/shraggy)
* [Will Sun](https://github.com/muirheadmaster)

## Acknowledgements
We would like to thank [James Wootton](https://github.com/quantumjim) for valuable suggestions and feedback. Our code closely follows his `RepetitionCode` structure in [`qiskit.ignis.verification.topological_codes`](https://qiskit.org/documentation/apidoc/verification.html#topological-codes), and his tutorials closely guided our initial explorations.

We'd also like to thank [Doug McClure](https://github.com/dtmcclure) for advising us on helpful details of the IBM hardware.

## References
- [Surface Codes: Towards Practical Large-Scale Quantum Computation](https://arxiv.org/abs/1208.0928)
- [Stabilizer Codes and Quantum Error Correction](https://arxiv.org/pdf/quant-ph/9705052.pdf)
- [Multi-path Summation for Decoding 2D Topological Codes](https://quantum-journal.org/wp-content/uploads/2018/10/q-2018-10-19-102.pdf)
- [Qiskit Textbook - Introduction to Quantum Error Correction using Repetition Codes](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html)
