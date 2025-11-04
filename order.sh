#!/bin/bash

START=0.1
INCREMENT=0.1
END=1.60
RUNS=5

echo "Temp,Order"
for TEMP in $(seq $START $INCREMENT $END)
do

  for (( i=1; i<=$RUNS; i++ ))
  do
    
    ORDER_VALUE=$(python LebwohlLasher_numba.py 1000 20 $TEMP 0 1 | awk -F 'Order: ' '{print $2}' | awk '{print $1}' | tr -d ',')

    echo "$TEMP,$ORDER_VALUE"
  done
  
done