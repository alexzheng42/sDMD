#!/bin/bash

pTime=(5000 25000 45000 65000 85000 105000 125000 145000 165000 185000 205000)
pTemp=(0.90  0.85  0.80  0.75  0.70   0.65   0.60   0.55   0.50   0.45   0.40) 

echo "Now processing ..."
head -3 ./input/parameter.txt
echo ""

./sDMD -f ala10.gro -i ./input -o ./

sed -i -e "1s/new/continue/" ./input/parameter.txt

for n in 0 1 2 3 4 5 6 7 8 9
	do
		sed -i -e "2s/${pTime[$n]}/${pTime[$n+1]}/" ./input/parameter.txt
		sed -i -e "3s/${pTemp[$n]}/${pTemp[$n+1]}/" ./input/parameter.txt

		echo "Now processing ..."
		head -3 ./input/parameter.txt
		echo ""

		./sDMD -f ala10.gro -i ./input -o ./
	done
