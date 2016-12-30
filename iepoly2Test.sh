#!/bin/bash
#Build
cd src
./buildIepoly2.sh
sleep 3
cp ./a.out ../
cd ..

#Execute
date && ./a.out && cat results/outputIepoly2.txt && date


#Generate graphs
cd src
/usr/bin/python3 iepoly2_graph.py

