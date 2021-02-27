# Qiskit Surface Code Encoder/Decoder
***Please check out our [presentation](https://docs.google.com/presentation/d/1HC5tQkvOcfl5lPDWy-l8ZvsWW3AGf6jUCmMu_sz-cDY/edit?usp=sharing) of this package to the Qiskit Advocate community at the November 2020 Qiskit Advocate Demo Session.***


---

# Setup

## 0. Getting setup with git version control and running bash scripts.

It may be worthwhile for Windows users to install a [github for windows](https://desktop.github.com/), which should come with Git Shell. 

If you are a Mac user, you *might* need XCode's [Command Line Package](https://developer.apple.com/library/archive/technotes/tn2339/_index.html) which allows for command line development in OSX. 

*Note: you do not need the XCode app to get these developer tools.* 

In order to install the package, run `xcode-select --install` in your terminal. If you get an error, you can troubleshoot and learn more about the package [here](https://developer.apple.com/opensource/).

## 1. Cloning the Repo

Run `git clone https://github.com/yaleqc/qiskit_surface_codes.git` to clone this repo. 

## 2. Install Miniconda (optional, but highly recommended)

As with many projects, it may be useful to set up a package manager **and** environment manager. [Miniconda](https://docs.conda.io/en/latest/miniconda.html#:~:text=Miniconda%20is%20a%20free%20minimal,zlib%20and%20a%20few%20others.) is a free, minimal installer for [conda](https://docs.conda.io/en/latest/) which servers as **both** a package and environment manager! 

You don't have to use miniconda, but we highly recommend it. Follow these [instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) for your operating system. Then, you should be able to run `conda -h` in your bash shell successfully. If this doesn't work, make sure to add miniconda to your [path](https://developers.google.com/earth-engine/guides/python_install-conda#windows_4).

## 3. Set up conda environment

Once you are able to successfully run `conda -h`, you can `cd`  (change directory) into the directory into which you just cloned the qiskit surface code project. 
Then, run `conda env create --file dependencies/requirements.yml` to install all necessary dependencies and create an environment called `qiskit-surface-code-env` using which you can test and develop the qiskit surface codes package!

*Note: The alternative is for you to install dependencies (python packages) manually until you are able to smoothly run the tutorials in this repo.*

## 4. Run the tutorials

Next, just activate your environment using `conda activate qiskit-surface-code-env`, `cd` into the project repo directory, and run `jupyter notebook`. 

Then, you should be able to open up the tutorial notebooks and run them without issues.


---

# About / History
## IBM Qiskit - Summer Jam Hackathon 2020

Quantum computation is an inherently noisy process. Scalable quantum computers will require fault-tolerance to implement useful computation. There are many proposed approaches to this, but one promising candidate is the family of *topological quantum error correcting codes*.

Currently, the [`qiskit.ignis.verification.topological_codes`](https://qiskit.org/documentation/apidoc/verification.html#topological-codes) module provides a general framework for QEC and implements one specific example, the *repetition code*.

For the hackathon, our team [Erwin's Tigers](#team) implemented a **surface code encoder and decoder** for Qiskit Ignis. We hope that this implementation will be useful to other Qiskitters and will inspire others to continue building out the `topological_codes` module into a diverse family.

## Usage
Inspired by the [Qiskit Textbook](https://qiskit.org/textbook/ch-quantum-hardware/error-correction-repetition-code.html), we've written a full set of [jupyter notebook tutorials](https://github.com/Phionx/qiskit_surface_codes/tree/master/tutorials), which are the best way to get up to speed. They detail both the API and the gritty implementation details -- please check them out!

## Background

<p align="center">
<img width="300" alt="surface code teaser" src="tutorials/img/error_set.jpg">
</p>

Surface codes are a type of CSS code, consisting of pairwise commuting X and Z stabilizers made of Pauli gates. It defines a logical state on a 2 by 2 lattice made of quantum bits with the stabilizers X<sub>1</sub> X<sub>2</sub> Z<sub>1</sub> Z<sub>2</sub>.

The code is based on the earlier theoretical idea of a *toric code*, with periodic boundary conditions instead of open boundary conditions. This has been shown to be largely identical, but embedding a surface code on an actual device is much easier.

## Implementation

In general, we try to follow the existing structure of [`qiskit.ignis.verification.topological_codes`](https://qiskit.org/documentation/apidoc/verification.html#topological-codes). The code is implemented separately here but is able to easily be merged into Ignis.

There are two main interfaces — corresponding to the encoder and decoder, respectively:

### `SurfaceCode` in [`topological_codes.circuits`](topological_codes/circuits/xxzz.py)

`SurfaceCode(d, T)` generates a `QuantumCircuit` for creating a logical state and measuring stabilizers. The class is parameterized with the code distance `d` (which should be odd) and the number of syndrome measurement rounds `T` (usually `T = d`). This class also handles parsing of the physical device readout into a form suitable for decoding. Please see the [encoder tutorial](tutorials/1_surface_code_encoding.ipynb) for a full walkthrough.
<p align="center">
<img width="615" alt="circuit" src="https://user-images.githubusercontent.com/293681/86277098-23ed3400-bba4-11ea-8305-c6d19eb73899.png">
</p>

### `GraphDecoder` in [`topological_codes.fitters`](topological_codes/fitters.py/xxzz.py)

`GraphDecoder(d, T)` implements minimum-weight perfect matching (MWPM) on the syndrome measurements of the physical circuit. This class is similar to the existing `GraphDecoder` for repetition codes, but introduces a new framework to handle the 2D lattice.

Parsed readout from the device is used to generate graphs of *error chains* in time and space, which decode syndrome measurements into the most likely sequence of qubit flips over time. Please see the [decoder tutorial](tutorials/2_surface_code_decoding.ipynb) for a full walkthrough.

<p align="center">
<img width="615" alt="matching graph" src="https://user-images.githubusercontent.com/293681/86277350-8ba37f00-bba4-11ea-9560-02d5ea3167cd.png">
</p>

## Future Directions
The scope of the project is quite large, so we focused on completing a "minimum viable product" during the hackathon. However, there are many areas which we'd like to explore going forward. A few immediate ones:

* Expand `SurfaceCode(d, T).circuits` a full set of logical states (1, +, -) -- and ultimately logical gate operations for computation.
* Full benchmark of the physical-error to logical-error probabilities to determine the error correction threshold.
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
