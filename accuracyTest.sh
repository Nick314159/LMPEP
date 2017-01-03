#!/bin/bash
#Build
cd src
./buildAccuracy.sh
sleep 3
cp ./a.out ../
cd ..

#Execute
date && ./a.out && cat results/outputAccuracy.csv && date


#Generate graphs
cd src
py=`which python`
$py accuracy_table.py

#Open graphs
cd ..
gnome-open results/accuracy_table.png
