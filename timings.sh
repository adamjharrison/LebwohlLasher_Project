#!/bin/bash
RUNS=5

echo "Size,Time"
for PARAM in 10 25 50 75 100 250 500 750 1000 2500 5000 7500 10000
do

  for (( i=1; i<=$RUNS; i++ ))
  do
    
    TIME=$(python LebwohlLasher.py 1000 $PARAM 0.65 0 1 | awk -F 'Time: ' '{print $2}' | awk '{print $1}' | tr -d ',')

    echo "$PARAM,$TIME"
  done
  
done