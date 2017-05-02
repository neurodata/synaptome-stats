#!/usr/bin/env bash
python getF0.py Ex12R76 data/Ex12R76/resultVol_1.csv resultVol_1 &
python getF0.py Ex12R76 data/Ex12R76/resultVol_2.csv resultVol_2 &
wait
python getF0.py Ex12R76 data/Ex12R76/resultVol_3.csv resultVol_3 &
python getF0.py Ex12R76 data/Ex12R76/resultVol_6.csv resultVol_6 &
python getF0.py Ex12R76 data/Ex12R76/resultVol_7.csv resultVol_7 &
wait
python getF0.py Ex12R76 data/Ex12R76/resultVol_8.csv resultVol_8 &
python getF0.py Ex12R76 data/Ex12R76/resultVol_9.csv resultVol_9 &
