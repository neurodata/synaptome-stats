#!/usr/bin/env bash
python getF0.py Ex3R43C1 data/Ex3R43C1/resultVol_1.csv resultVol_1 &
python getF0.py Ex3R43C1 data/Ex3R43C1/resultVol_2.csv resultVol_2 &
python getF0.py Ex3R43C1 data/Ex3R43C1/resultVol_3.csv resultVol_3 &
wait
python getF0.py Ex3R43C1 data/Ex3R43C1/resultVol_6.csv resultVol_6 &
python getF0.py Ex3R43C1 data/Ex3R43C1/resultVol_7.csv resultVol_7 &
python getF0.py Ex3R43C1 data/Ex3R43C1/resultVol_8.csv resultVol_8 &
python getF0.py Ex3R43C1 data/Ex3R43C1/resultVol_9.csv resultVol_9 &
