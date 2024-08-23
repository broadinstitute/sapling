Sapling
=======

Simple simulations of pandemic-scale sequencing data with known parameters.

Overview
--------

Given a population model, a sampling strategy, a substitution model and a genome size, Sapling simulates
the sequencing data that would be seen for a sample of N individuals under a coalescent model with no selection.
First, we pick the sampling times of all N individuals.  Then, using the population model, we run a coalescent
simulation to generate a genealogy for them.  Then, using the substitution model and genome size, we generate
a random sequence for the most-recent common ancestor (MRCA).  Finally, we evolve the sequence along the
tree using the Gillespie algorithm (similar to [phastSim](https://github.com/NicolaDM/phastSim)).

The output consists of:
1) an info file with input model parameters
2) a MAPLE file for the sequences
3) a Newick tree annotated with a reference sequence at the root and with mutation events along the branches.

At the moment, there is no site rate heterogeneity and no sequences are missing any data.

Sapling was written specifically to easily generate large simulated datasets to benchmark
[Delphy](https://github.com/broadinstitute/delphy).  It is similar, but more limited than,
[phastSim](https://github.com/NicolaDM/phastSim) and [Seq-Gen](https://github.com/rambaut/Seq-Gen).  The main
differences are that both ancestry and sequence evolution are simulated, not just the latter.  Like phastSim, but unlike
Seq-Gen, the full set of sequences is never materialized (internally, sapling build up an Explicit Mutation-Annotated
Tree), so very large datasets can be generated with modest resources.  By outputting the results in MAPLE format, the
output can also be modest in size, provided that the simulation parameters correspond to a typical genomic epidemiology
dataset (i.e., very few mutations between a sequence and its closest common ancestor with the remaining sequences).

System Requirements
-------------------
Sapling is written in generic C++20 and can be compiled as a native, standalone command-line program.  Build instructions can be found in [`INSTALL.md`](INSTALL.md).  Sapling was developed and tested under Linux (Ubuntu 22.04.4 LTS, x86-64).

Credits and Acknowledgements
----------------------------

Sapling is developed in the [Sabeti Lab](https://www.sabetilab.org/) at the [Broad
Institute](https://www.broadinstitute.org/).

Sapling draws a lot of inspiration from:

- [Delphy](https://github.com/broadinstitute/delphy)
- [phastSim](https://github.com/NicolaDM/phastSim)
- [Seq-Gen](https://github.com/rambaut/Seq-Gen)

Copyright (c) 2022-2024 Broad Institute, Inc.  See [LICENSE](LICENSE) for details.
