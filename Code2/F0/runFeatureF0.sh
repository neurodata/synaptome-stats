#!/usr/bin/env bash
TOKEN1=$1 #"kristina15"
INPUT=$2 #"testLocations.csv"
OUTPUT1=$3 #"outJ.h5"
BF=$4 #"5" ## Buffer size around each center point

echo $TOKEN1
echo $INPUT
echo $OUTPUT1
echo $BF

##getF*y <token> <input> <output> <buffer> 
getF0.py $TOKEN1 $INPUT $OUTPUT1 $BF
