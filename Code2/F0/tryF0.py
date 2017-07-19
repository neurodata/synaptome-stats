import sys
import os
from multiprocessing import Pool
import csv
import datetime
import numpy as np
from numpy import genfromtxt
import h5py
import ndio
import ndio.remote.neurodata as neurodata

def chan(TOKEN):
    nd = neurodata()
    ### Pull channel names
    channels = nd.get_channels(TOKEN)
    
    ## Filter channels corresponding to images only
    chan = []
    for key in channels:
        if channels[key]["channel_type"] == "image" or channels[key]["channel_type"] == "oldchannel":
            chan.append(key)
    
    columnNames = [str(key) for key in chan]
    return(columnNames)


def F0(L, bf, TOKEN, channels):
    nd = neurodata()
    cols = []
    xLocal = L[0]
    yLocal = L[1]
    zLocal = L[2]
    print(xLocal, yLocal, zLocal)
    syn = []
    cube = []
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
            f0 = np.sum(np.asarray(gg))
        except: 
            syn.append(-1000)
            r = bf*2 +1
            out = np.empty([r,r,r])
            out.fill(-1000)
        else:
            syn.append(f0)
            cube.append(gg)
            
    outCubes = np.asarray(cube)
    F0dat = np.array(syn, dtype = 'int32')
    
    dout = {"cubes": outCubes, "F0": F0dat}
    return(F0dat)

def F0dict(L, bf, TOKEN, channels):
    nd = neurodata()
    cols = []
    xLocal = L[0]
    yLocal = L[1]
    zLocal = L[2]
    print(xLocal, yLocal, zLocal)
    syn = {}
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
            f0 = np.sum(np.asarray(gg))
        except: 
            syn.update({ch : -1000})
        else:
            syn.update({ch : f0})

    #F0dat = np.array(syn, dtype = 'int32')

    return(F0dat)


def F0k15(l):
    ## This is a specific channel order to make downstream
    ## analysis easier.
    channels = ["synapsinR_7thA", "synapsinGP_5thA", "VGluT1_3rdA", 
         "VGluT1_8thA", "VGluT2_2ndA", "psd_8thA", "GluR2_2ndA",
         "NMDAR1_6thA", "NR2B_9thA", "NOS_9thA", "Synpod_3rdA",
         "GAD_6thA", "VGAT_5thA", "PV_1stA", "gephyrin_1stA",
         "GABARa1_4thA", "GABABR1_3rdA", "VGluT3_1stA", "CR1_2ndA",
         "TH_5thA", "VAChT_4thA", "tubulin_8thA", "dapi_1stA"]

    return(F0(l, 5, "kristina15", channels))

def F0weiler(l):
    ## This is a specific channel order to make downstream
    ## analysis easier.

    channels = ['Synapsin1_2', 'vGluT1_3', 'vGluT2_2', 'PSD95_1',
            'GluR1_4', 'GluR2_8', 'GluR4_8', 'NR2A_1', 'NR2B_3',
            'Synaptopodin_6', 'GABAARa1_7', 'GAD2_3', 'Gephyrin_2',
            'PV25_5', 'vGAT_4', 'vGAT_5', 'vGAT_6', 'GFP_4', 'YFP_1',
            'Arc_5', 'Calbindin_7', 'DAPI_1', 'DAPI_2', 'DAPI_3',
            'DAPI_4', 'DAPI_5', 'DAPI_6', 'DAPI_7', 'DAPI_8' ]

    return(F0(l, 5, "Ex10R55", channels))


def main(n, TOKEN, INPUT, OUTPUT, BF):
    nd = neurodata()
    ## This is a specific channel order to make downstream
    ## analysis easier.
    
    loc = genfromtxt(INPUT, delimiter=',', skip_header = 0).tolist()
    L = np.around(np.asarray(loc), decimals = 0).astype(int)

    #print(L)
    #out = F0k15(L[0])
    #print(out)
    with Pool(n) as p:
        out = p.map(F0weiler,L)

    h5f0OUT = h5py.File(OUTPUT+".h5", 'w')
    h5f0OUT.create_dataset(TOKEN + "_F0", data = np.transpose(out))
    h5f0OUT.create_dataset("Locations", data = np.transpose(L))
    h5f0OUT.close()

    channels = ['Synapsin1_2', 'vGluT1_3', 'vGluT2_2', 'PSD95_1',
            'GluR1_4', 'GluR2_8', 'GluR4_8', 'NR2A_1', 'NR2B_3',
            'Synaptopodin_6', 'GABAARa1_7', 'GAD2_3', 'Gephyrin_2',
            'PV25_5', 'vGAT_4', 'vGAT_5', 'vGAT_6', 'GFP_4', 'YFP_1',
            'Arc_5', 'Calbindin_7', 'DAPI_1', 'DAPI_2', 'DAPI_3',
            'DAPI_4', 'DAPI_5', 'DAPI_6', 'DAPI_7', 'DAPI_8' ]

    np.savetxt(OUTPUT+".csv", out, delimiter=",", fmt = "%d", header =
            ','.join(channels), comments = "")
    np.savetxt("testLocationsOUT.csv",L, delimiter=",")
    #return(out)


if __name__ == '__main__':
    n = 4
    ## input parameters
    TOKEN = sys.argv[1] ## e.g. "Ex10R55"  "kristina15"
    INPUT = sys.argv[2] ## a csv file of locations in (x, y, z) order.
    OUTPUT = sys.argv[3] ## A base name for output
    BF     = sys.argv[4] ## Buffer size around each center point

    nd = neurodata()
    main(n, TOKEN, INPUT, OUTPUT, BF)


