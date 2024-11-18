#!/bin/bash

echo "Execution Times (n=2 to n=100)"
for ((k=2; k<=100; k++))
do
    sed -i "9s/.*/#define n $k/" code1.c
    bash make.sh

done

echo "Execution timing completed and stored in out.dat"

