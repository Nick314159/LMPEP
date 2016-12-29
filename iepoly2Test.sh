#!/bin/bash
cd src
./buildIepoly2.sh
sleep 3
cp ./a.out ../
cd ..

date && ./a.out && cat results/outputIepoly2.txt && date



