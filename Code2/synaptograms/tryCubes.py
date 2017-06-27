#!/usr/bin/python3
import sys
import os
import itertools
from functools import partial
from multiprocessing import Pool
from multiprocessing import cpu_count
import csv
import datetime
import numpy as np
from numpy import genfromtxt
import h5py
import ndio
import ndio.remote.neurodata as neurodata

nd = neurodata("synaptomes.neurodata.io")

def chan(TOKEN):
    nd = neurodata("synaptomes.neurodata.io")
    #nd = neurodata()
    ### Pull channel names
    channels = nd.get_channels(TOKEN)
    
    ## Filter channels corresponding to images only
    chan = []
    for key in channels:
        if channels[key]["channel_type"] == "image":
            chan.append(key)
    
    columnNames = [str(key) for key in chan]
    return(columnNames)

def links(token, loc, channels=None, col=None):
    ## I know there's probably a much better way to do this. 
    if channels == None:
        channels = chan(token)

    if col == None:
        col = [1 for i in range(len(channels))]

    x = str(loc[0])
    y = str(loc[1])
    z = str(loc[2])

    bSt = "https://viz-dev.boss.neurodata.io/#!{'layers':{"
    bEnd = "}_'navigation':{'pose':{'position':{'voxelSize':[100_100_70]_'voxelCoordinates':["
    bEnd = bEnd + x + "_" + y +"_" + z + "]}}_'zoomFactor':70}}"
    
    chMid = "':{'type':'image'_'source':'boss://https://api.boss.neurodata.io/" + token + "/image/"
    chEnd = lambda col :"?window=0,10000'_'opacity':0.42_'color':" + str(col) + "}"
    
    chs = ["'" + ch + chMid + ch + chEnd(cl) for ch, cl in zip(channels, col)]
        
    sep = "_"
    sep = sep.join(chs)

    out = bSt + sep + bEnd
    return(out)


#def cube(L, bf, TOKEN, channels):
def cube(d):
    L = d['l'] 
    bf = int(d['bf'])
    TOKEN =d['token'] 
    channels = d['channels']

    cols = []
    xLocal = L[0]
    yLocal = L[1]
    zLocal = L[2]
    print(xLocal, yLocal, zLocal)
    cubeCh = []
    for ch in channels:
        try:
            queryGlobal = {
                    'token': TOKEN,
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
        except: 
            r = bf*2 +1
            out = np.empty([r,r,r])
            out.fill(-1000)
        else:
            cubeCh.append(gg)
            out = np.asarray(cubeCh)

    return(out)


def cubeK15(d):

    channels=["synapsinR_7thA", "synapsinGP_5thA", "VGluT1_3rdA",
            "VGluT1_8thA", "VGluT2_2ndA", "psd_8thA", "GluR2_2ndA",
            "NMDAR1_6thA", "NR2B_9thA", "NOS_9thA", "Synpod_3rdA",
            "GAD_6thA", "VGAT_5thA", "PV_1stA", "gephyrin_1stA",
            "GABARa1_4thA", "GABABR1_3rdA", "VGluT3_1stA", "CR1_2ndA",
            "TH_5thA", "VAChT_4thA", "tubulin_8thA", "dapi_1stA"]

    ## The default channel argument is a specific 
    ## channel order to make downstream analysis easier.

    di = {
         "l": d['l'], 
         "bf":d['bf'], 
         "token":"kristina15",
         "channels":channels
         }

    ret = cube(di)
    return(ret)


def main(n, TOKEN, INPUT, BF):

    loc = genfromtxt(INPUT, delimiter=',', skip_header = 1).tolist()
    L = np.around(
          np.asarray(loc), decimals = 0).astype(np.int32)[...,0:3]

    di = []
    for a in range(len(L)):
        di.append({'l' : L[a], 'bf': BF})

    with Pool(n) as p:
        out = p.map(cubeK15, di)

    arOut = {
            "cubes": np.asarray(out, dtype = np.int32), 
            "loc"  : L
            }

    #return(np.asarray(out))
    return(arOut)



def test():
    n = 5
    ## input parameters
    TOKEN = "kristina15"
    INPUT = "c1121.csv"
    OUTPUT = "zooOUTPUT"
    BF = 5

    nd = neurodata("synaptomes.neurodata.io")
    n = 5

    ## Get cubes and return as np.array
    channels = ["synapsinR_7thA", "synapsinGP_5thA", "VGluT1_3rdA", 
         "VGluT1_8thA", "VGluT2_2ndA", "psd_8thA", "GluR2_2ndA",
         "NMDAR1_6thA", "NR2B_9thA", "NOS_9thA", "Synpod_3rdA",
         "GAD_6thA", "VGAT_5thA", "PV_1stA", "gephyrin_1stA",
         "GABARa1_4thA", "GABABR1_3rdA", "VGluT3_1stA", "CR1_2ndA",
         "TH_5thA", "VAChT_4thA", "tubulin_8thA", "dapi_1stA"]

    out = main(n, TOKEN, INPUT, BF)

    ## Save data
    h5f0OUT = h5py.File(OUTPUT+".h5", 'w')

    h5f0OUT.create_dataset(TOKEN + "_cubes", data = np.transpose(out["cubes"]))
    h5f0OUT.create_dataset("Channels", data = np.string_(channels))
    h5f0OUT.create_dataset("Locations", data = np.transpose(out["loc"]))
    h5f0OUT.close()

    return(out)

def chans():
    channels = ["synapsinR_7thA", "synapsinGP_5thA", "VGluT1_3rdA", 
         "VGluT1_8thA", "VGluT2_2ndA", "psd_8thA", "GluR2_2ndA",
         "NMDAR1_6thA", "NR2B_9thA", "NOS_9thA", "Synpod_3rdA",
         "GAD_6thA", "VGAT_5thA", "PV_1stA", "gephyrin_1stA",
         "GABARa1_4thA", "GABABR1_3rdA", "VGluT3_1stA", "CR1_2ndA",
         "TH_5thA", "VAChT_4thA", "tubulin_8thA", "dapi_1stA"]
    return(channels)


if __name__ == '__main__':

    TOKEN = sys.argv[1] ## e.g. "Ex10R55"  "kristina15"
    INPUT = sys.argv[2] ## a csv file of locations (x, y, z) order.
    BF    = sys.argv[3] ## Buffer size around each center point
    OUTPUT= sys.argv[4] ## A base name for output

    nd = neurodata("synaptomes.neurodata.io")
    #nd = neurodata() 

    ## Number of processes to use
    #n = multiprocessing.cpu_count()
    n = 5

    ## Get cubes and return as np.array
    channels = ["synapsinR_7thA", "synapsinGP_5thA", "VGluT1_3rdA", 
         "VGluT1_8thA", "VGluT2_2ndA", "psd_8thA", "GluR2_2ndA",
         "NMDAR1_6thA", "NR2B_9thA", "NOS_9thA", "Synpod_3rdA",
         "GAD_6thA", "VGAT_5thA", "PV_1stA", "gephyrin_1stA",
         "GABARa1_4thA", "GABABR1_3rdA", "VGluT3_1stA", "CR1_2ndA",
         "TH_5thA", "VAChT_4thA", "tubulin_8thA", "dapi_1stA"]

    out = main(n, TOKEN, INPUT, BF)

    ## Save data
    h5f0OUT = h5py.File(OUTPUT+".h5", 'w')

    h5f0OUT.create_dataset(TOKEN + "_cubes", data = np.transpose(out["cubes"]))
    h5f0OUT.create_dataset("Channels", data = np.string_(channels))
    h5f0OUT.create_dataset("Locations", data = np.transpose(out["loc"]))
    h5f0OUT.close()




