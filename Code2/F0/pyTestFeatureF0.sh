#!/usr/bin/env bash
TOKEN1="kristina15"
INPUT="locationTest5.csv"
OUTPUT1="foo999.h5"
BF="5" ## Buffer size around each center point

echo $TOKEN1
echo $INPUT
echo $OUTPUT1
echo $BF

##python tryF*y <token> <input> <output> <buffer> 
python3 tryF0.py $TOKEN1 $INPUT $OUTPUT1 $BF
