#!/usr/bin/env bash
TOKEN1=$1 #"kristina15"
INPUT=$2 #"testLocations.csv"
OUTPUT1=$3 #"outJ.h5"
OUTPUT2=$4 #"k15Synaptograms1000"
BF=$5 #"5" ## Buffer size around each center point

mkdir ~/data/outputs

echo $1 
echo $2
echo $3 
echo $4 
echo $5 

##python getCub*.py <token> <input> <output> <buffer> 
cd /home/meda/data

getCubes.py $TOKEN1 $INPUT outputs/$OUTPUT1 $BF
#
##Rscript k15Synaptogr*.R <input>  <token> <output> <buffer>
k15Synaptograms.R outputs/$OUTPUT1 $TOKEN1 outputs/$OUTPUT2 $BF

