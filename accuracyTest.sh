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
            echo 'Runs the accuracy test. Use -p to print result data. Use -d  to build with debug flags. Use -o to open graphs and table upon generation. -h for this help message'
            exit
            ;;
        --debug | -d)
            FLAGS=$FLAGS' -d'
            ;;
        --open | -o)
            OPEN=true
            ;;
         --print | -p)
            PRINT=true
            ;;
        *)
            PARAMATERS=$PARAMATERS" $1" 
    esac
    shift
done

#Build
cd src
./buildAccuracy.sh $FLAGS
sleep 3
cd ..

#Execute
echo "Starting Accuracy test at `date`" 
bin/accuracy.out $PARAMTERS
if $PRINT 
then 
cat results/outputAccuracy.csv
fi
echo "Finished Accuracy test at `date`" 

#Generate graphs
cd src
py=`which python`
$py accuracy_table.py

#Open graphs
if $OPEN 
then 
cd ..
xdg-open results/accuracy_table.png
fi
