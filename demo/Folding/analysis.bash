#!/bin/bash
# Analyze the folding demo: energy, hydrogen bonds, radius of gyration and
# RMSD to an ideal alpha-helix and beta-hairpin. -sum joins all 11 runs.
# Run from this folder after building sAnalysis in ../../analysis.
# Note: -ref file names are relative to the -path folder.

set -e
../../analysis/sAnalysis -path output/ -trj out_trj.gro -cnt out_cnt.txt -sum \
    -rPBC -HB -En -RG -RMSD -ref ../analysis/ref_alpha.gro ../analysis/ref_beta.gro
mv output/{Energy,HBInfo,RG,RMSD}.txt analysis/

# Plot the results (needs python3 with numpy and matplotlib).
(cd analysis && python3 plot.py) || echo "Plotting skipped: needs python3 with numpy and matplotlib."
