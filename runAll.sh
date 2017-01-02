#!/bin/bash
clear
echo 'Cleaning results dir for new run...'
./clean.sh
echo ' ------------------------------------------'
echo 'Running spoly...'
./spolyTest.sh 10 1000
echo ' ------------------------------------------'
echo 'Running gepoly...'
./gepolyTest.sh 10 100 10 100
echo ' ------------------------------------------'
echo 'Running iepoly1...'
./iepoly1Test.sh 10 100 10 100
echo ' ------------------------------------------'
echo 'Running iepoly2...'
./iepoly2Test.sh
echo ' ------------------------------------------'
echo 'Running accuracy...'
./accuracyTest.sh
