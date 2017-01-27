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
            echo 'Runs the tridiagonal test. Use -p to print result data. Use -d  to build with debug flags. Use -o to open graphs and table upon generation. -h for this help message'
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
cd src
./buildTri.sh $FLAGS
sleep 3
cd ..

#Execute
echo "Starting Tridiagonal test at `date`" 
./bin/tri.out $PARAMATERS 
#> /dev/null 2>&1
if $PRINT 
then
  cat results/outputTri.csv
fi
echo "Finished Tridiagonal test at `date`" 

#Generate graphs
cd src
py=`which python`
$py tri_table.py

#Open graphs
if $OPEN 
then 
cd ..
for file in `ls results/tri_*.pdf`
do
xdg-open $file
done
fi



