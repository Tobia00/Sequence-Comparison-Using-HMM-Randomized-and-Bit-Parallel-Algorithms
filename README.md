# Sequence-Comparison-Using-HMM-Randomized-and-Bit-Parallel-Algorithms

Biological Sequence Alignment Package
=====================================

This package contains a compact implementation of three approaches for biological
sequence comparison: randomized alignment (Monte Carlo and Las Vegas), Pair-HMM
alignment using the Viterbi algorithm, and a bit-parallel inspired global
alignment score comparison. The package is intended as a small user-manual style
submission for the computational biology algorithms project.

Installation
------------
No external Python libraries are required. Use Python 3.9 or newer.

1. Unzip the package.
2. Open a terminal in the unzipped folder.
3. Run the demo with:

   python SRC/demo.py

Usage
-----
The main scripts are located in the SRC directory:

- randomized_alignment.py: Monte Carlo and Las Vegas sequence alignment methods.
- pair_hmm_alignment.py: Pair-HMM alignment with Match, Insertion, and Deletion states.
- bit_parallel_alignment.py: Classical DP score and BitPAl-style score comparison.
- demo.py: Runs a short demonstration using the sample FASTA files in SRC/data.

To use your own FASTA files, import the functions from the SRC scripts or edit
SRC/demo.py and replace the demo FASTA paths with your own sequence file paths.

Demo
----
The demo uses two small FASTA files provided in SRC/data. It prints randomized
alignment results, Pair-HMM alignment metrics, and a DP/BitPAl-style score check.
This keeps the package portable and easy to test on any machine.
