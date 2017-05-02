#!/usr/bin/env bash
python getF0.py Ex12R75 data/Ex12R75/resultVol_1.csv resultVol_1 &
python getF0.py Ex12R75 data/Ex12R75/resultVol_2.csv resultVol_2 &
wait
python getF0.py Ex12R75 data/Ex12R75/resultVol_3.csv resultVol_3 &
python getF0.py Ex12R75 data/Ex12R75/resultVol_6.csv resultVol_6 &
wait
python getF0.py Ex12R75 data/Ex12R75/resultVol_7.csv resultVol_7 &
python getF0.py Ex12R75 data/Ex12R75/resultVol_8.csv resultVol_8 &
wait
python getF0.py Ex12R75 data/Ex12R75/resultVol_9.csv resultVol_9 &
