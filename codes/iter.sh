#!/bin/bash

echo "Execution Times (n=1 to n=1000)"
for ((k=200; k<=1000; k = k+50))
do
    sed -i "9s/.*/#define n $k/" code1.c
    bash make.sh

done

echo "Execution timing completed and stored in out.dat"

