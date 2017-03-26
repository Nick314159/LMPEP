#!/bin/bash

#Setup
FLAGS=''
PARAMETERS=''
while [ ! $# -eq 0 ]
do
    case "$1" in
        --help | -h)
            echo 'Runs the whole test. Use -p to print result data. Use -d  to build with debug flags. Use -o to open graphs and table upon generation. -h for this help message'
            exit
            ;;
        --debug | -d)
            FLAGS=$FLAGS' -d'
            ;;
         --print | -p)
            FLAGS=$FLAGS' -p'
            ;;
        --open | -o)
            FLAGS=$FLAGS' -o'
            ;;
    esac
    shift
done


clear
./bin/clean.sh
echo ' ------------------------------------------'
echo 'Running spoly...'
./spolyTest.sh $FLAGS 50 10000 5
echo ' ------------------------------------------'
echo 'Running gepoly...'
./gepolyTest.sh $FLAGS 20 400 50 2000 5
echo ' ------------------------------------------'
echo 'Running cc...'
./ccTest.sh $FLAGS 20 640 50 2000 5
echo ' ------------------------------------------'
echo 'Running iepoly1...'
./iepoly1Test.sh $FLAGS 20 400 50 2000 5
echo ' ------------------------------------------'
echo 'Running iepoly2...'
./iepoly2Test.sh $FLAGS 
echo ' ------------------------------------------'
echo 'Running accuracy...'
./accuracyTest.sh $FLAGS 
echo ' ------------------------------------------'
echo 'Running Tridiagonal...'
./triTest.sh $FLAGS 
echo ' ------------------------------------------'
echo 'Running Sample Tridiagonal...'
./sampleTest.sh $FLAGS 
