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
cd src
./buildIepoly2.sh $FLAGS
sleep 3
cd ..

#Execute
echo "Starting Iepoly2 test at `date`" 
bin/iepoly2.out $PARAMTERS
if $PRINT 
then 
cat results/outputIepoly2-*.txt
fi
echo "Finished Iepoly2 test at `date`" 



#Generate graphs
cd src
py=`which python`
$py iepoly2_graph.py

#Open graphs
if $OPEN 
then 
cd ..
gnome-open results/iepoly2*.pdf
fi
