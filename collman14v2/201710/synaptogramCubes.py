#!/usr/bin/env python3
###
###  
###
### Jesse Leigh Patsolic 
### 2017 <jpatsol1@jhu.edu>
### S.D.G 
#
import argparse
import math
from intern.remote.boss import BossRemote
from intern.resource.boss.resource import *
import configparser
#import grequests # for async requests, conflicts with requests somehow
import requests
import numpy as np
from numpy import genfromtxt
import h5py
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import shutil
import blosc
from IPython.core.debugger import set_trace
import sys
import os
import itertools
from functools import partial
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
import csv
import datetime
import morton
import toolbox



def main(COLL_NAME, EXP_NAME, COORD_FRAME, LOCATIONS, BF = [5,5,5],
         CHAN_NAMES=None, num_threads = 6, CONFIG_FILE= None):

    L = genfromtxt(LOCATIONS,  delimiter=',', skip_header = 0, dtype='int32').tolist()
    m = morton.Morton(dimensions=3, bits=64)
    Lzo = sorted([m.pack(L[i][0], L[i][1], L[i][2]) for i in range(len(L))])

    loc = np.asarray([m.unpack(l) for l in Lzo])

    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    TOKEN = config['Default']['token']
    boss_url = ''.join( ( config['Default']['protocol'],'://',config['Default']['host'],'/v1/' ) )
    
    #intern
    rem = BossRemote(CONFIG_FILE)

    cfr = rem.get_project(CoordinateFrameResource(COORD_FRAME))
    
    ext = {'x': [cfr.x_start, cfr.x_stop], 'y': [cfr.y_start, cfr.y_stop], 'z': [cfr.z_start, cfr.z_stop]}

    ## remove coordinates that are outside buffer area. 
    x_bnd = [np.all([x[0]-BF[0] >= ext['x'][0] , x[0]+BF[0] < ext['x'][1]]) for x in loc] 
    y_bnd = [np.all([y[1]-BF[1] >= ext['y'][0] , y[1]+BF[1] < ext['y'][1]]) for y in loc] 
    z_bnd = [np.all([z[2]-BF[2] >= ext['z'][0] , z[2]+BF[2] < ext['z'][1]]) for z in loc] 

    xyzbnd = np.asarray([np.all([x_bnd[i], y_bnd[i], z_bnd[i]]) for i in range(len(x_bnd))])

    ## Select only locations that are buffered away from the edge.
    loc = loc[xyzbnd] 

    ## Compute the buffered ranges 
    x_rng = [[x[0]-BF[0], x[0]+BF[0]+1] for x in loc] 
    y_rng = [[y[1]-BF[1], y[1]+BF[1]+1] for y in loc] 
    z_rng = [[z[2]-BF[2], z[2]+BF[2]+1] for z in loc]

    ChanList = []
    for ch in CHAN_NAMES:
        di = [{
              'rem': rem,
              'ch_rsc':
                ChannelResource(ch,COLL_NAME,EXP_NAME,'image',datatype='uint8'),
              'ch'  : ch,
              'res' : 0,
              'loc' : loc[i],
              'xrng': x_rng[i],  
              'yrng': y_rng[i],  
              'zrng': z_rng[i],  
              'bf'  : BF 
             } for i in range(len(loc))]
   
        with ThreadPool(num_threads) as tp:
            out = tp.map(toolbox.getCube, di)
            print(ch) ##DEBUG
            sys.stdout.flush() #DEBUG

        ChanList.append(np.asarray(out))

    outArray = np.asarray(ChanList)
    ## outArray will be a numpy array in [channel, synapse, z, y, x] order
    ## this is C-ordering
    return(outArray, loc) 
## END MAIN

def driveMain():
    COLL_NAME      = 'collman' 
    EXP_NAME       = 'collman14v2' 
    COORD_FRAME    = 'collman_collman14v2'
    LOCATIONS_FILE = 'meda_plots_tight/hmcC5synaptogramLocations.csv'
    BF             = [108,108,5]
    OUTPUT         = 'blankTest'
    CONFIG_FILE    = 'config.ini'

    #CHAN_NAMES = ['DAPI1st', 'DAPI2nd', 'DAPI3rd', 'GABA488',
    #            'GAD647', 'gephyrin594', 'GS594', 'MBP488', 'NR1594',
    #            'PSD95_488', 'Synapsin647', 'VGluT1_647']
    #CHAN_NAMES = ['synapsin', 'PSDr']
    CHAN_NAMES = ['EM10K', 'bIIItubulin', 'DAPI_2nd', 'DAPI_3rd',
		  'GABA', 'GAD2', 'VGAT', 'gephyrin', 'NR1', 'PSDr', 'synapsin', 'VGluT1']

    cubes, locs = main(COLL_NAME, EXP_NAME, COORD_FRAME, LOCATIONS_FILE, BF = BF,
         CHAN_NAMES=CHAN_NAMES, num_threads = 6, CONFIG_FILE= 'config.ini')

    import pdb; pdb.set_trace()
    f0 = toolbox.F0(cubes)
    #toolbox.mainOUT(f0, CHAN_NAMES, OUTPUT)
    toolbox.toh5(EXP_NAME, OUTPUT + '.h5', CHAN_NAMES, cubes, locs, F0 = f0)

    return(cubes, locs)
## END driveMain


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 
        'Download synapse cubes using the BOSS')
    parser.add_argument('-C', help='Valid collection id',
            type = str, metavar='C', default='weiler14')
    parser.add_argument('-E', help='Valid experiment id',
            type = str, metavar='E', default='Ex2R18C1')
    parser.add_argument('-F', help='valid coordinate frame', 
            type = str, metavar='F', default='weiler14_Ex2R18C1')
    parser.add_argument('-B', help='integer xyz list for buffer around center',
            nargs = 3, type = int, metavar='B', default= [5,5,5])
    parser.add_argument('-L', help='csv file of locations '
            'in xyz order', type = str, metavar='L', required=True)
    parser.add_argument('-O', help='output filename',
            type = str, metavar='O', required=True,
            default = 'output')
    parser.add_argument('--con', help='user config file for BOSS'
            'authentication', type = str, metavar='con', required=True)
    
    args = parser.parse_args()

    COLL_NAME      = args.C
    EXP_NAME       = args.E
    COORD_FRAME    = args.F
    LOCATIONS_FILE = args.L 
    BF             = args.B
    OUTPUT         = args.O
    CONFIG_FILE    = args.con

    #rem = BossRemote(CONFIG_FILE)
    #CHAN_NAMES = rem.list_channels(COLL_NAME, EXP_NAME)

    CHAN_NAMES = ['bIIItubulin', 'DAPI_2nd', 'DAPI_3rd', 
                  'GABA', 'GAD2', 'VGAT', 'gephyrin', 
                  'NR1', 'PSDr', 'synapsin', 'VGluT1'] 

    cubes, locs = main(COLL_NAME, EXP_NAME, COORD_FRAME, 
                       LOCATIONS_FILE, BF = BF, 
                       CHAN_NAMES=CHAN_NAMES, 
                       num_threads = 6, CONFIG_FILE= CONFIG_FILE)

    f0 = toolbox.F0(cubes)

    toolbox.mainOUT(f0, CHAN_NAMES, OUTPUT)
    toolbox.toh5(EXP_NAME, OUTPUT + '.h5', CHAN_NAMES, cubes, locs, F0 = f0)
    print('Done!')

