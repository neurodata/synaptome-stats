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
    nd = neurodata("synaptomes.neurodata.io")
    ### Pull channel names
    channels = nd.get_channels(TOKEN)
    
    ## Filter channels corresponding to images only
    chan = []
    for key in channels:
        if channels[key]["channel_type"] == "image":
            chan.append(key)
    
    columnNames = [str(key) for key in chan]
    return(columnNames)


def F0(L, bf, TOKEN, channels):
    nd = neurodata("synaptomes.neurodata.io")
    cols = []
    xLocal = L[0]
    yLocal = L[1]
    zLocal = L[2]
    print(xLocal, yLocal, zLocal)
    syn = []
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
        else:
            syn.append(f0)

    F0dat = np.array(syn, dtype = 'int32')
    return(F0dat)


def F0k15(l):
    channels = ['CR1_2ndA', 'GABABR1_3rdA', 'GABARa1_4thA', 
        'GAD_6thA', 'GluR2_2ndA', 'NMDAR1_6thA', 
        'NOS_9thA', 'NR2B_9thA', 'PV_1stA',
        'Synpod_3rdA', 'TH_5thA', 'VAChT_4thA', 'VGAT_5thA',
        'VGluT1_3rdA', 'VGluT1_8thA', 'VGluT2_2ndA', 'VGluT3_1stA',
        'dapi_1stA', 'gephyrin_1stA', 'psd_8thA', 'synapsinGP_5thA',
        'synapsinR_7thA', 'tubulin_8thA']
    return(F0(l, 5, "kristina15", channels))


def main():
    channels = ['CR1_2ndA', 'GABABR1_3rdA', 'GABARa1_4thA', 
        'GAD_6thA', 'GluR2_2ndA', 'NMDAR1_6thA', 
        'NOS_9thA', 'NR2B_9thA', 'PV_1stA',
        'Synpod_3rdA', 'TH_5thA', 'VAChT_4thA', 'VGAT_5thA',
        'VGluT1_3rdA', 'VGluT1_8thA', 'VGluT2_2ndA', 'VGluT3_1stA',
        'dapi_1stA', 'gephyrin_1stA', 'psd_8thA', 'synapsinGP_5thA',
        'synapsinR_7thA', 'tubulin_8thA']
    ## input parameters
    #TOKEN = sys.argv[1] ## e.g. "Ex10R55"  "kristina15"
    #INPUT = sys.argv[2] ## a csv file of locations in (x, y, z) order.
    #OUTPUT = sys.argv[3] ## A base name for output
    #BF     = sys.argv[4] ## Buffer size around each center point
    TOKEN = 'kristina15'
    INPUT = './locationTest48bad.csv'
    OUTPUT = './testOUT'
    BF = 5

    loc = genfromtxt(INPUT, delimiter=',', skip_header = 1).tolist()
    L = np.around(np.asarray(loc), decimals = 0).astype(int)

    #print(L)
    #out = F0k15(L[0])
    #print(out)
    with Pool(8) as p:
        out = p.map(F0k15,L)

    h5f0OUT = h5py.File(OUTPUT+".h5", 'w')
    h5f0OUT.create_dataset(TOKEN + "_F0", data = np.transpose(out))
    h5f0OUT.create_dataset("Locations", data = np.transpose(L))
    h5f0OUT.close()

    np.savetxt(OUTPUT+".csv", out, delimiter=",")
    np.savetxt("testLocationsOUT.csv",L, delimiter=",")
    #return(out)


if __name__ == '__main__':
    main()


