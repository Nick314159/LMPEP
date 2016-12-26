#!/bin/bash
cd src
./build.sh
sleep 3 
cp ./a.out ../

cd ..
date && ./a.out $1 $2 $3 && cat results/output.csv && date

cd src/AMVW/doubleshift
make
cd tests
./test $2
cd ../../../../

