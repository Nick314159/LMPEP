#!/bin/bash
cd src
./build.sh
sleep 3
cp ./a.out ../
cd ..

cd src/AMVW/doubleshift
make > /dev/null 2>&1
cd tests
date 
./test $2 > /dev/null 2>&1
cd ../../../../
cat results/outputAMVW.csv
date

date && ./a.out $1 $2 $3 && cat results/output.csv && date



