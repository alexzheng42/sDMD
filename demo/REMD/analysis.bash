#!/bin/bash
# Analyze the REMD demo: per-replica energy, hydrogen bonds, RMSD to an ideal
# alpha-helix and beta-hairpin and the potential-energy map, then WHAM with the
# RMSD as reaction coordinate (heat capacity, free energy, combined energy map).
# Run from this folder after building sAnalysis in ../../analysis.
# Note: -ref file names are relative to the -path folder.

set -e
../../analysis/sAnalysis -path output/ -trj out_trj.gro -cnt out_cnt.txt \
    -rPBC -HB -En -RMSD -PotMap -ref ../analysis/ref_alpha.gro ../analysis/ref_beta.gro \
    -REMD -reNo 12 -rmsd -rate 50 -bin 100
mv output/{Energy,HBInfo,RMSD,PotMap}*.txt output/{PMFMap_T,FreeEnergy,Probability,CvPotT,CombinedPotMap}*.txt analysis/ 2>/dev/null || true

# Plot the results (needs python3 with numpy and matplotlib).
(cd analysis && python3 plot.py) || echo "Plotting skipped: needs python3 with numpy and matplotlib."
