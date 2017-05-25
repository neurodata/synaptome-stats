#!/usr/bin/env python
import csv
import datetime
import ndio
import numpy as np
import h5py
import scipy.linalg
from numpy import genfromtxt
import ndio.remote.neurodata as neurodata
import os
import sys

### Load user defined functions
execfile("/bin/functions.py")

## input parameters
TOKEN1 = sys.argv[1] ## e.g. "Ex10R55"  "kristina15"
INPUT  = sys.argv[2] ## a csv file of locations in (x, y, z) order.
OUTPUT = sys.argv[3] ## A base name for output
BF     = sys.argv[4] ## Buffer size around each center point

## TESTING PARAMETERS
# TOKEN1 = "kristina15" #REMOVE
# INPUT  = "testLocations.csv"#REMOVE
# OUTPUT = "outJ" #REMOVE
# BF     = 5 #REMOVE

bf = int(BF)
nd = neurodata(hostname = "synaptomes.neurodata.io")
output = OUTPUT


###  Read in location data from input
loc = genfromtxt(INPUT, delimiter=',', skip_header = 1).tolist()

### Calculate global dimensions of the volume
size = nd.get_image_size(TOKEN1)
globalMaxs = [size[0], size[1], size[2]]

### Pull channel names
channels = nd.get_channels(TOKEN1)

## Filter channels corresponding to images only
chan = []
for key in channels:
    if channels[key]["channel_type"] == "image":
        chan.append(key)

channels = chan

columnNames = [str(key) for key in channels]


### Filter locations that are far enough away from the
### boundary of the volume that a full cube may be obtained
L1 = []
for i in range(len(loc)):
    if inRange(loc[i], globalMaxs, bf):
        L1.append(loc[i])
        
L = np.around(np.asarray(L1), decimals = 0).astype(int)


### For each location grab all channels and compute sum statistic
cols = []
for i in range(len(L)):
    xLocal = L[i][0]
    yLocal = L[i][1]
    zLocal = L[i][2]
    syn = []
    for ch in chan:
        queryGlobal = {
                'token': TOKEN1,
                'channel': ch,
                'x_start': xLocal - bf,
                'x_stop': xLocal + bf +1,
                'y_start': yLocal - bf,
                'y_stop': yLocal + bf +1,
                'z_start': zLocal - bf,
                'z_stop': zLocal + bf +1,
                'resolution': 0
                }
        gg = nd.get_cutout(**queryGlobal)
        f0 = np.sum(np.asarray(gg))
        syn.append(f0)
    cols.append(syn)

F0dat = np.array(cols)


### Write output to hdf5 file
h5f0OUT = h5py.File(output, 'w')
h5f0OUT.create_dataset(TOKEN1 + "_F0", data = np.transpose(F0dat))
h5f0OUT.create_dataset("Locations", data = L)
h5f0OUT.create_dataset("Channels", data = columnNames)
h5f0OUT.close()
