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
cd bin
./buildIepoly1.sh $FLAGS
sleep 3
cd ..

#Execute
echo "Starting Iepoly1 test at `date`" 
bin/iepoly1.out $PARAMATERS > /dev/null 2>&1 
if $PRINT 
then 
cat results/outputIepoly1*.csv
fi
echo "Finished Iepoly1 test at `date`" 

#Generate graphs
cd python
py=`which python`
$py iepoly1_graph.py

#Open graphs
if $OPEN 
then 
cd ..
for file in `ls results/iepoly_times_*.pdf`
do
xdg-open $file
done
fi
