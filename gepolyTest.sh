#!/bin/bash
cd src
./buildGepoly.sh
sleep 3
cp ./a.out ../
cd ..

date && ./a.out $1 $2 $3 $4 && cat results/outputGepoly.csv && date



