#!/bin/bash
#Build
cd src
./buildIepoly1.sh
sleep 3
cp ./a.out ../
cd ..

#Execute
date && ./a.out $1 $2 $3 $4 > /dev/null 2>&1 && cat results/outputIepolySize1.csv && cat results/outputIepolyDegree1.csv&& date

#Generate graphs
cd src
py=`which python`
$py iepoly1_graph.py
