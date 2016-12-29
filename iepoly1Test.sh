#!/bin/bash
cd src
./buildIepoly1.sh
sleep 3
cp ./a.out ../
cd ..

date && ./a.out $1 $2 $3 $4 && cat results/outputIepolySize1.csv && cat results/outputIepolyDegree1.csv&& date



