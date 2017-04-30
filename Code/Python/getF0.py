
# coding: utf-8

# # resultsVol_1: Getting data
# 
# Anish ran synapse detection on the chessboard dataset and has provided locations of potential synapses.
# Here, we pull data from those locations and export it into a nice dataset for further analysis.

# In[8]:

import csv
import ndio ### MUST HAVE VERSION 1.1.5
import numpy as np
import h5py
import scipy.linalg
from numpy import genfromtxt
import ndio.remote.neurodata as neurodata
import os
import sys

nd = neurodata()
local = False

TOKEN1 = (sys.argv[1])##"Ex10R55"
INPUT = (sys.argv[2])
NAME = (sys.argv[3])

bf = 5

output = TOKEN1 + "_" + NAME + ".h5"


# ### User defined functions

# In[9]:

execfile("functions.py")


# In[10]:

public_tokens = nd.get_public_tokens()


# In[11]:

TOKEN1 in public_tokens # Should *definitely* be true


# In[12]:

size = nd.get_image_size(TOKEN1)
channels = nd.get_channels(TOKEN1)

xGlobal = [0,size[0]]
yGlobal = [0,size[1]]
zGlobal = [0,size[2]]

globalMaxs = [size[0], size[1], size[2]]


# In[13]:

chan = []
for key in channels:
    chan.append(key)


# ## Import resultsVol_1
# Import the location data from "resultVol_1.csv" making sure that each location is not too 
# close to the boundary of the data and round to nearest integer.
# And that the points are far enough away from the boundary.

# In[14]:

loc = genfromtxt(INPUT, delimiter=',').tolist()

L1 = []
for i in range(len(loc)):
    if inRange(loc[i], globalMaxs, bf):
        L1.append(loc[i])

L = np.around(np.asarray(L1), decimals = 0).astype(int)


# ## Extract 11x11x11 cube around each location 
# ### 11x11x11 cube is used for features F0-F3, with the buffer needed for features F4-F5
# then cast to vector of length 11^3 and output as an array dim = (n,1331)
# 
# #### N.B. This does take a while

# In[18]:

cols = []
for ch in chan:
    queryGlobal = {
            'token': TOKEN1,
            'channel': ch,
            'x_start': xGlobal[0],
            'x_stop': xGlobal[1],
            'y_start': yGlobal[0],
            'y_stop': yGlobal[1],
            'z_start': zGlobal[0],
            'z_stop': zGlobal[1],
            'resolution': 0
            }
    gg = nd.get_cutout(**queryGlobal)
    #print "got",ch
    
    OUT = []
    
    for i in range(len(L)):
        cube11 = getCube(gg, bf, L[i])
        F = f0(cube11)
        OUT.append(F)
        
    out = np.asarray(OUT)
    cols.append(out)
    #print "got F0 for",ch
    


# In[48]:

columnDat = np.transpose(np.asarray(cols))
columnNames = np.string_(chan)


# ## Save as HDF5

# In[49]:

h5fOUT = h5py.File(output, 'w')
h5fOUT.create_dataset(TOKEN1, data=columnDat)
h5fOUT.create_dataset('Locations', data=np.transpose(L))
h5fOUT.close()

print "saved " + output


