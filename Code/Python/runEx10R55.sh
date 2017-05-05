#!/usr/bin/env bash
python getF0.py Ex10R55 LocationDataXYZ/Ex10R55/resultVol_1xyz.csv resultVol_1 &
wait
python getF0.py Ex10R55 LocationDataXYZ/Ex10R55/resultVol_2xyz.csv resultVol_2 &
wait
python getF0.py Ex10R55 LocationDataXYZ/Ex10R55/resultVol_3xyz.csv resultVol_3 &
python getF0.py Ex10R55 LocationDataXYZ/Ex10R55/resultVol_6xyz.csv resultVol_6 &
wait
python getF0.py Ex10R55 LocationDataXYZ/Ex10R55/resultVol_7xyz.csv resultVol_7 &
python getF0.py Ex10R55 LocationDataXYZ/Ex10R55/resultVol_8xyz.csv resultVol_8 &
wait
python getF0.py Ex10R55 LocationDataXYZ/Ex10R55/resultVol_9xyz.csv resultVol_9 &
