#!/bin/bash
RUNS=5

echo "Size,Time"
for PARAM in 10	18	32	56	100	178	316	562	1000	1778	3162	5623	10000
do

  for (( i=1; i<=$RUNS; i++ ))
  do
    
    TIME=$(python LebwohlLasher_mpi4py.py 1000 $PARAM 0.65 0 1 | awk -F 'Time: ' '{print $2}' | awk '{print $1}' | tr -d ',')

    echo "$PARAM,$TIME"
  done
  
done