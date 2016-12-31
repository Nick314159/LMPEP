#!/bin/bash
#Configure
cd results
resultsDir=`pwd`
cd ..

#Build
cd src
./buildSpoly.sh
sleep 3
cp ./a.out ../
cd AMVW/doubleshift
./configure.sh $resultsDir
make > /dev/null 2>&1

#Execute
cd tests
date 
./spoly $1 $2 $3 > /dev/null 2>&1
cd ../../../../
cat results/outputAMVW.csv
date
echo
date && ./a.out $1 $2 $3 && cat results/output.csv && date

#Generate graphs
cd src
/usr/bin/python spoly_graph.py



