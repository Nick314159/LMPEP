#!/bin/bash
#Build
cd src
./buildIepoly2.sh
sleep 3
cp ./a.out ../
cd ..

#Execute
date && ./a.out && cat results/outputIepoly2-*.txt && date


#Generate graphs
cd src
py=`which python`
$py iepoly2_graph.py

#Open graphs
cd ..
gnome-open results/iepoly2*.pdf
