#!/usr/bin/env bash
TOKEN1=$1 #"kristina15"
INPUT=$2 #"testLocations.csv"
OUTPUT1=$3 #"outJ.h5"
OUTPUT2=$4 #"k15Synaptograms1000"
BF=$5 #"5" ## Buffer size around each center point
MEANS=$6
SDS=$7

echo $1 
echo $2
echo $3 
echo $4 
echo $5 
echo $6
echo $7

##python getCub*.py <token> <input> <output> <buffer> 
cd /home/meda/data

tryCubes.py $TOKEN1 $INPUT $BF $OUTPUT1 
#
##Rscript k15Synaptogr*.R <input>  <token> <output> <buffer>
k15Synaptograms.R $OUTPUT1.h5 $TOKEN1 $OUTPUT2 $BF $MEANS $SDS


