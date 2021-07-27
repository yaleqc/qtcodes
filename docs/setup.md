# Setup

## 0. Getting setup with git version control and running bash scripts.

It may be worthwhile for Windows users to install a [github for windows](https://desktop.github.com/), which should come with Git Shell.

If you are a Mac user, you *might* need XCode's [Command Line Package](https://developer.apple.com/library/archive/technotes/tn2339/_index.html) which allows for command line development in OSX.

*Note: you do not need the XCode app to get these developer tools.*

In order to install the package, run `xcode-select --install` in your terminal. If you get an error, you can troubleshoot and learn more about the package [here](https://developer.apple.com/opensource/).

## 1. Cloning the Repo

Run `git clone https://github.com/yaleqc/qiskit_topological_codes.git` to clone this repo.

## 2. Install Miniconda (optional, but highly recommended)

As with many projects, it may be useful to set up a package manager **and** environment manager. [Miniconda](https://docs.conda.io/en/latest/miniconda.html#:~:text=Miniconda%20is%20a%20free%20minimal,zlib%20and%20a%20few%20others.) is a free, minimal installer for [conda](https://docs.conda.io/en/latest/) which servers as **both** a package and environment manager!

You don't have to use miniconda, but we highly recommend it. Follow these [instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html#regular-installation) for your operating system. Then, you should be able to run `conda -h` in your bash shell successfully. If this doesn't work, make sure to add miniconda to your [path](https://developers.google.com/earth-engine/guides/python_install-conda#windows_4).

## 3. Set up conda environment

Once you are able to successfully run `conda -h`, you can `cd`  (change directory) into the directory into which you just cloned the qiskit topological code project.
Then, run `conda env create --file dependencies/requirements.yml` to install all necessary dependencies and create an environment called `qiskit-topological-codes-env` using which you can test and develop the qiskit topological codes package!

*Note: The alternative is for you to install dependencies (python packages) manually until you are able to smoothly run the tutorials in this repo.*

## 4. Run the tutorials

Next, just activate your environment using `conda activate qiskit-topological-codes-env`, `cd` into the project repo directory, and run `jupyter lab`.

Then, you should be able to open up the tutorial notebooks and run them without issues.