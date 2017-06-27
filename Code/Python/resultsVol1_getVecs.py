
# coding: utf-8

# # resultsVol_1: Getting data
#
# Anish ran synapse detection on the chessboard dataset and has provided locations of potential synapses.
# Here, we pull data from those locations and export it into a nice dataset for further analysis.

# In[1]:

import csv
import ndio
import numpy as np
import h5py
import scipy.linalg
from numpy import genfromtxt
import ndio.remote.neurodata as neurodata

nd = neurodata()


# ### User defined functions

# In[3]:

def buff( x, bf ):
    return [x - bf, x + bf]

inRange = lambda x,y,bf: bf <= x[0] <= (y[0]-1-bf) and bf <= x[1] <= (y[1]-1-bf) and bf <= x[2] <= (y[2]-1-bf)

def get11Bool( loc , big, bf = 5 ):
    xyz = np.around(loc, decimals=0).astype(int)
    x = [xyz[0] - bf,xyz[0] + bf + 1]
    y = [xyz[1] - bf,xyz[1] + bf + 1]
    z = [xyz[2] - bf,xyz[2] + bf + 1]
    #
    Bool = min(x) >= 0 and min(y) >= 0 and min(z) >=0 and    max(x) < big.shape[0] and max(y) < big.shape[1] and max(z) < big.shape[2]
    #
    if Bool:
        cube = big[(xyz[0] - bf):(xyz[0] + bf + 1):1,                   (xyz[1] - bf):(xyz[1] + bf + 1):1,                   (xyz[2] - bf):(xyz[2] + bf + 1):1]
        return cube

def get11( loc , big, bf = 5 ):
    xyz = np.around(loc, decimals=0).astype(int)
    #
    cube = big[(xyz[0] - bf):(xyz[0] + bf + 1):1,             (xyz[1] - bf):(xyz[1] + bf + 1):1,             (xyz[2] - bf):(xyz[2] + bf + 1):1]
    return cube

def cube2vec( cube ):
    vec = np.concatenate(np.concatenate(cube))
    return(vec)


# In[3]:

public_tokens = nd.get_public_tokens()

TOKEN1 = "Ex10R55"
#CHANNEL = ["vGAT_6","Synapsin1_2"]
CHANNEL = "Synapsin1_2"

TOKEN1 in public_tokens # Should *definitely* be true


# In[4]:

size = nd.get_image_size(TOKEN1)
xGlobal = [0,size[0]]
yGlobal = [0,size[1]]
zGlobal = [0,size[2]]

globalMaxs = [size[0], size[1], size[2]]
bf = 3


# In[5]:

queryGlobal = {
        'token': TOKEN1,
        'channel': CHANNEL,
        'x_start': xGlobal[0],
        'x_stop': xGlobal[1],
        'y_start': yGlobal[0],
        'y_stop': yGlobal[1],
        'z_start': zGlobal[0],
        'z_stop': zGlobal[1],
        'resolution': 0
        }

gg = nd.get_cutout(**queryGlobal)


# ## Import resultsVol_1
# Import the location data from "resultVol_1.csv" making sure that each location is not too
# close to the boundary of the data and round to nearest integer.

# In[ ]:

loc = genfromtxt('../../data/Anish_csvfiles/Ex10R55/resultVol_1.csv', delimiter=',').tolist()

L1 = []
for i in range(len(loc)):
    if inRange(loc[i], globalMaxs, 5):
        L1.append(loc[i])

L = np.around(np.asarray(L1), decimals = 0).astype(int)


# ## Save data to an HDF5 file for local use.

# In[ ]:

#h5f = h5py.File('Ex10R55_Synapsin1_2.h5', 'w')
#h5f.create_dataset('Ex10R55_Synapsin1_2', data=gg)
#h5f.close()


# ## Extract 11x11x11 cube around each location
# then cast to vector of length 11^3 and output as an array dim = (n,1331)

# In[9]:

OUT = []
for i in range(len(L)):
    OUT.append(cube2vec(get11(L[i],gg, bf)))

out = np.asarray(OUT)


# In[12]:

L


# ## Save data to HDF5 and csv

# In[16]:

h5fOUT = h5py.File('7cubeEx10R55_Synapsin1_2_vecs.h5', 'w')
h5fOUT.create_dataset('7cEx10R55_Synapsin1_2_resultVol1', data=np.transpose(out))
h5fOUT.create_dataset('Locations', data=np.transpose(L))
h5fOUT.close()


# In[17]:

np.savetxt("Ex10R55_Synapsin1_2_Locations.csv", L)
np.savetxt("Ex10R55_Synapsin1_2_vecs.csv", out)


# In[ ]:


F0 = []
for i in range(len(OUT)):
    F0.append(f0(OUT[i]))


F1 = []
for i in range(len(OUT)):
    F1.append(f2(OUT[i]))


A = np.array([np.sqrt((i-(bf+1))**2 + (j-(bf+1))**2 + (k-(bf+1))**2) for i in range(1,2*bf+2) for j in range(1,2*bf+2) for k in range(1,2*bf+2)])

B = A.reshape((7,7,7))



