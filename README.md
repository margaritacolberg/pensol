# PenSol: Protein Folding in a Penetrating Solvent

[![build workflow status](https://github.com/margaritacolberg/pensol/actions/workflows/build.yml/badge.svg)](https://github.com/margaritacolberg/pensol/actions/workflows/build.yml?query=branch:main)
[![format workflow status](https://github.com/margaritacolberg/pensol/actions/workflows/format.yml/badge.svg)](https://github.com/margaritacolberg/pensol/actions/workflows/format.yml?query=branch:main)

PenSol simulates a protein folding in a penetrating solvent. The solvent
particles are not explicitly defined. Instead, only the random forces that
occur on the surface of the protein to mimic solvent-protein collisions are
simulated. Ten sample transitions of the protein crambin are chosen between
different pairs of intermediate states, and the transition rate (`krate`) and
the population of the unbonded state at equilibrium (`Pu`) are calculated for
each.

For more information, see the [preprint](https://arxiv.org/abs/2310.13223).

## Directories

See top level comment in each file for more details as to what each file does.
Some top level comments also contain examples of how to run each file, if
complicated. The pathways given in these examples assume user is currently in
the `pensol/crambin` directory. Download `crambin.zip` from
[Releases](https://github.com/margaritacolberg/pensol/releases).

  * `src`: C++ source code for PenSol program

  * `examples`: JSON inputs which specify protein configuration for running
    PenSol

  * `tools`: Python scripts for running PenSol and extracting data such as
    the transition rate

  * `crambin`: HybridMC JSON and HDF5 inputs for running specific transitions
    for crambin between intermediate states; dir inside which all simulations
    are carried out

## C++ Program Details

### Prerequisite Software

  * C++17 compiler with support for C++20
    [bit library](https://en.cppreference.com/w/cpp/header/bit)

    see bit operations under C++20 [compiler
    support](https://en.cppreference.com/w/cpp/compiler_support/20)

    tested with GCC 9.3.0 on CentOS 7.9

  * Minimum CMake version is 3.13.1

  * [Ninja](https://ninja-build.org/)

  * HDF5 version 1.10 or later

  * Boost

### Commands for Compiling

To compile PenSol for development and testing:

```
make debug
```

To perform unit test with verbose output:

```
cd debug
ctest -V
```

To compile PenSol for running simulations:

```
make release
```

## Python Program Details

### Create JSON Inputs

The JSON file needed to run a transition between two intermediate states of
crambin uses inputs from a master JSON file for crambin inside `examples`
called `crambin.json`, and the HybridMC JSON file specific to the transition in
`crambin`. The inputs from these two files must be combined into one JSON file
for each transition, which can be done using `update_json.py`. This step must
be completed before the program is run.

### Run Program

To run a simulation for one transition on a computing cluster, use
`run_diff_eps.py`. To run a simulation for one transition on a local computer,
use `run.py`.

### Calculate PenSol `krate` and `Pu`

To calculate the `krate` and `Pu` for the PenSol model, use a non-zero epsilon
value such as 3.0 when submitting a simulation using `run_diff_eps.py`. When
all transitions at this epsilon are completed, use `run_krate_solvent.py` to
get `krate` and `Pu`. This will create a csv file in the home directory called
`krate_solvent.csv`, which stores `krate` and `Pu`. Alternatively, to obtain
`krate` and `Pu` for a single transition, enter transition specific folder, and
run `krate_solvent.py`.

### Calculate HybridMC Entropy and Mean First Passage Times (MFPT)

To get the `krate` for the HybridMC model, the biased entropy and MFPT must be
calculated beforehand. This can be done in each `crambin/hybridmc_*` directory.
To calculate biased entropy, use `get_error_s_bias.py` for one or more `nboot`
(for crambin, 5 `nboot` is sufficient). To calculate MPFT, use `mfpt.py`. For
crambin, bootstraps are not needed in the MPFT calculation. If there are
several HDF5 files in `crambin/hybridmc_*` due to `nboot` of
`get_error_s_bias.py` being greater than 1, use the first HDF5 file to
calculate MFPT.

### Calculate HybridMC `krate` and `Pu`

To calculate the `krate` and `Pu` for the HybridMC model, set epsilon to zero
when submitting a simulation using `run_diff_eps.py`. When all transitions at
this epsilon are completed, use `run_vacf_or_msd.py` with `msd.py` to calculate
the mean squared displacement (MSD) for each transition. To get `krate` and
`Pu` for the HybridMC model for all transitions, use `run_krate_hybridmc.py`.
This will create a csv file in the main directory called `krate_hybridmc.csv`,
which stores `krate` and `Pu`. Alternatively, to obtain `krate` and `Pu` for a
single transition, enter transition specific folder and run `plot_msd.py`,
followed by `krate_hybridmc.py`.

If the system consists of a single Brownian particle, either `plot_msd.py` or
`plot_vacf.py` can be used to obtain the short time diffusion coefficient. In
the latter case, `run_vacf_or_msd.py` would need to be carried out with
`vacf.py` to calculate the velocity autocorrelation function (VACF).

### Rate Ratios

To obtain a ratio of hybridmc to solvent `krate` at different epsilons for a
single transition, run `rate_ratios.py`.

### Calculate the Diffusion Coefficient

To calculate the diffusion coefficient based on the values of the parameters
given in `crambin.json`, use `calc_diff.py`.
