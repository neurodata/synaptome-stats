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

#TOKEN1 = "kristina15" #REMOVE
#INPUT  = "testLocations.csv"#REMOVE
#OUTPUT = "outJ" #REMOVE
#BF     = 5 #REMOVE

TOKEN1 = sys.argv[1] ## e.g. "Ex10R55"  "kristina15"
INPUT  = sys.argv[2] ## a csv file of locations in (x, y, z) order.
OUTPUT = sys.argv[3] ## A base name for output
BF     = sys.argv[4] ## Buffer size around each center point

bf = int(BF)
nd = neurodata(hostname = "synaptomes.neurodata.io")
output = OUTPUT

### Load user defined functions
execfile("/bin/functions.py")


size = nd.get_image_size(TOKEN1)
globalMaxs = [size[0], size[1], size[2]]

channels = nd.get_channels(TOKEN1)

chan = []
for key in channels:
    if channels[key]["channel_type"] == "image":
        chan.append(key)


channels = chan

columnNames = [str(key) for key in channels]

loc = genfromtxt(INPUT, delimiter=',', skip_header = 1).tolist()


# In[4]:

L1 = []
for i in range(len(loc)):
    if inRange(loc[i], globalMaxs, bf):
        L1.append(loc[i])
        
L = np.around(np.asarray(L1), decimals = 0).astype(int)


# In[5]:

synCube = []
for i in range(len(L)):
    cubeCh = []
    for ch in channels:
        xLocal = L[i][0]
        yLocal = L[i][1]
        zLocal = L[i][2]
        queryGlobal = { 
            'token': TOKEN1,
            'channel': ch,
            'x_start': xLocal - bf,
            'x_stop': xLocal + bf +1,
            'y_start': yLocal - bf,
            'y_stop': yLocal + bf +1,
            'z_start': zLocal - bf,
            'z_stop': zLocal + bf + 1,
            'resolution': 0   
            }
        
        gg = nd.get_cutout(**queryGlobal)
        cubeCh.append(gg)
        out = np.asarray(cubeCh)
    synCube.append(out)


# In[6]:

### Save as an hdf5 file
h5fOUT = h5py.File(output, 'a')
h5fOUT.create_dataset(TOKEN1, data=np.transpose(synCube))
h5fOUT.create_dataset('Locations', data=np.transpose(L))
h5fOUT.create_dataset('Channels', data=columnNames)
h5fOUT.close()

