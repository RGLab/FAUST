#!/bin/bash

for i in `seq 1 3`;
do
    for j in `seq 1 2`;
    do
	for k in `seq 1 2`;
	do
	    for l in `seq 1 3`;
	    do
		/bin/bash ./launchSim.sh $i $j $k $l
	    done
	done
    done
done 
