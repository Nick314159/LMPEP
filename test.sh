#!/bin/bash
cd src
./build.sh
sleep 3 
cp ./a.out ../

cd ..
date && ./a.out $1 $2 $3 && cat results/output.csv && date




#date=`date +"%d-%m-%Y"`
#output="testn$1$date.log"
#for ((n = 1; n <= $1; n++)); do
#	degree="1^$n" | bc`
#	start=`date`
#	echo =========================>> $output
#	echo         DEGREE $degree>> $output
#	echo =========================>> $output
#	echo Start: $start >> $output 	
#	
#	./a.out 1 $degree 2
#	wait $!
#	
#	end=`date`
#	cat results/output.csv >> $output
#	echo End: $end >> $output 	
#done
