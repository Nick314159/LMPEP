#!/bin/bash
#Setup
FLAGS=''
OPEN=false
PRINT=false
PARAMETERS=''
while [ ! $# -eq 0 ]
do
    case "$1" in
        --help | -h)
            echo 'Runs the spoly test. Use -p to print result data. Use -d  to build with debug flags. Use -o to open graphs and table upon generation. -h for this help message'
            exit
            ;;
        --debug | -d)
            FLAGS=$FLAGS' -d'
            ;;
         --print | -p)
            PRINT=true
            ;;
        --open | -o)
            OPEN=true
            ;;
        *)
            PARAMATERS=$PARAMATERS" $1" 
    esac
    shift
done

#Configure
cd results
resultsDir=`pwd`
cd ..

#Build
cd bin
./buildSpoly.sh $FLAGS
sleep 3
cd ../src/AMVW/doubleshift
./configure.sh $resultsDir
make > /dev/null 2>&1

#Execute
cd tests
echo "Starting Spoly AMVW test at `date`" 
./spoly $PARAMATERS > /dev/null 2>&1
cd ../../../../
if $PRINT 
then
  cat results/outputAMVW.csv
fi
echo "Finished Spoly AMVW test at `date`" 

echo

echo "Starting Spoly test at `date`" 
bin/spoly.out $PARAMATERS
if $PRINT 
then 
  cat results/outputSpoly.csv
fi
echo "Finished Spoly test at `date`" 

#Generate graphs
cd python
py=`which python`
$py spoly_graph.py

#Open graphs
if $OPEN 
then 
cd ..
for file in `ls results/spoly_*.pdf`
do
xdg-open $file
done
fi



