#!/usr/bin/env bash
TOKEN1="kristina15"
INPUT="locationTest5.csv"
OUTPUT1="test.h5"
BF="5" ## Buffer size around each center point

echo $TOKEN1
echo $INPUT
echo $OUTPUT1
echo $BF

##python tryF*y <token> <input> <output> <buffer> 
python tryF0.py $TOKEN1 $INPUT $OUTPUT1 $BF
