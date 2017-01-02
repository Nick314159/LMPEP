#!/bin/bash
#Build
cd src
./buildAccuracy.sh
sleep 3
cp ./a.out ../
cd ..

#Execute
date && ./a.out && cat results/outputAccuracy.txt && date


#Generate graphs
cd src
py=`which python`
$py accuracy_table.py

