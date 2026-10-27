#!/bin/bash
# REMD demo: 12 replicas of ALA10 between 236.5 and 416 K, 200000 time units each.
# Run from this folder after "make all" in ../../src. Uses 12 CPU cores.

set -e
cp ../../src/sDMD ../../src/sREMD ../../src/sServer .
mkdir -p input output
cp parameter.txt ala10.gro REMDConfig.txt input/
[ -e input/Library_Ding ] || ln -s ../../../input/Library_Ding input/Library_Ding

./sREMD -f input/REMDConfig.txt -args -i input -o output -f ala10.gro -end
