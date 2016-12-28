#!/bin/bash
cd src
./buildSpoly.sh
sleep 3
cp ./a.out ../
cd ..

date && ./a.out $1 $2 $3 && cat results/outputGepoly.csv && date



