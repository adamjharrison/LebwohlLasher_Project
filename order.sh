#!/bin/bash

START=0.05
INCREMENT=0.05
END=1.60
RUNS=5

PARAM=${1:-"0"}
echo "Temp,Order"
for TEMP in $(seq $START $INCREMENT $END)
do

  for (( i=1; i<=$RUNS; i++ ))
  do
    
    ORDER_VALUE=$(python LebwohlLasher.py 1000 20 $TEMP $PARAM 1 | awk -F 'Order: ' '{print $2}' | awk '{print $1}' | tr -d ',')

    echo "$TEMP,$ORDER_VALUE"
  done
  
done