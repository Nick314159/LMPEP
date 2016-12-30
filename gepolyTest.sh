#!/bin/bash
#Build and setup
cd src
./buildGepoly.sh
sleep 3
cp ./a.out ../
cd ..

#Execute
date && ./a.out $1 $2 $3 $4 > /dev/null 2>&1 && cat results/outputGepolySize.csv && cat results/outputGepolyDegree.csv&& date

#Generate graphs
cd src
/usr/bin/python3 gepoly_graph.py


