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
./clean.sh
echo ' ------------------------------------------'
echo 'Running spoly...'
./spolyTest.sh $FLAGS 10 1000
echo ' ------------------------------------------'
echo 'Running gepoly...'
./gepolyTest.sh $FLAGS 10 100 10 100
echo ' ------------------------------------------'
echo 'Running iepoly1...'
./iepoly1Test.sh $FLAGS 10 100 10 100
echo ' ------------------------------------------'
echo 'Running iepoly2...'
./iepoly2Test.sh $FLAGS 
echo ' ------------------------------------------'
echo 'Running accuracy...'
./accuracyTest.sh $FLAGS 
