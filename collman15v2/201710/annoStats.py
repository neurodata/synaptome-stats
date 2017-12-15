#!/usr/bin/env python3
###
###  
###
### Jesse Leigh Patsolic 
### 2017 <jpatsol1@jhu.edu>
### S.D.G 
#
import argparse
from intern.remote.boss import BossRemote
from intern.resource.boss.resource import *
from intern.utils.parallel import block_compute
import configparser
import requests
import numpy as np
from numpy import genfromtxt
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
import toolbox


def main(COLL_NAME, EXP_NAME, COORD_FRAME,
         CHAN_NAMES=None, num_threads = 4, CONFIG_FILE= 'config.ini'):

    bf = [5,245,245] # in z,y,x order

    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    TOKEN = config['Default']['token']
    boss_url = ''.join( ( config['Default']['protocol'],'://',config['Default']['host'],'/v1/' ) )
    #print(boss_url)
    #'https://api.boss.neurodata.io/v1/'
    
    #intern
    rem = BossRemote(CONFIG_FILE)

    cf = CoordinateFrameResource(str(COLL_NAME + '_' + EXP_NAME))
    cfr = rem.get_project(cf)
    anno_res = ChannelResource('annotation', COLL_NAME, EXP_NAME, 'annotation', datatype='uint64')
    
    ex = {'x': cfr.x_stop, 'y': cfr.y_stop, 'z': cfr.z_stop}
    
    blocks = block_compute(0,ex['x'],0,ex['y'],0,ex['z'],
               origin = (0,0,0), block_size = (512, 512, 16))

    rid = []
    for b in blocks:
      rid = rid + rem.get_ids_in_region(anno_res, 0, b[0], b[1], b[2], [0,1])

    u = np.unique(np.asarray(rid))

    ## bounding box for annotation_i
    bb = [rem.get_bounding_box(anno_res, 0,ui, 'tight') for ui in u]

    for i in range(len(bb)):
        bb[i]["id"] = u[i]

    A = [(rem.get_cutout(
          anno_res, 0, bb[i]["x_range"], 
          bb[i]["y_range"], bb[i]["z_range"], 
          id_list = [bb[i]['id']]) != 0).astype(int) 
        for i in range(len(bb))] 

    #Bmeans = [np.int32(np.round(np.mean(np.asarray(np.where(A[i] == True)),1))) for i in range(len(A))]
    Bmeans = [np.int32(np.round(np.mean(np.asarray(np.where(A[i] == 1)),1))) for i in range(len(A))]
    
    Bglobal = []
    for i in range(len(bb)):
        ad1 = np.asarray([bb[i]['z_range'][0], bb[i]['y_range'][0], bb[i]['x_range'][0]])
        Bglobal.append(Bmeans[i] + ad1)
        
    ColMin = np.asarray(bf)
    ColMax = np.asarray([ex['z'] - (bf[0] + 1),  # The z index is inclusive
                         ex['y'] - (bf[1] + 1), 
                         ex['x'] - (bf[2] + 1)])

    m = [Bglobal[i] >= ColMin for i in range(len(Bglobal))]
    M = [Bglobal[i] <= ColMax for i in range(len(Bglobal))]
    mm = [np.all(m[i]) for i in range(len(m)) ]
    MM = [np.all(M[i]) for i in range(len(M)) ]

    Bcon = []
    con = [np.asarray(mm[j] and MM[j]) for j in range(len(mm))]
    for i in range(len(Bglobal)):
        if con[i]: 
            Bcon.append(Bglobal[i])

       
    if CHAN_NAMES is None:
                CHAN_NAMES = ['DAPI1st', 'DAPI2nd', 'DAPI3rd',
        	      'GABA488', 'GAD647', 'gephyrin594', 'GS594', 'MBP488',
                      'NR1594', 'PSD95_488', 'Synapsin647', 'VGluT1_647']
                #CHAN_NAMES = ['bIIItubulin', 'DAPI_2nd', 'DAPI_3rd',
                #        'GABA', 'GAD2', 'gephyrin', 'NR1', 'PSDr',
                #        'synapsin', 'VGAT', 'VGluT1'] 

    ChanList = []

    ## for getting masked 
    #for ch in CHAN_NAMES:
    #    di = [{
    #          'rem': rem,
    #          'ch_rsc':
    #            ChannelResource(ch,COLL_NAME,EXP_NAME,'image',datatype='uint8'),
    #          'ch'  : ch,
    #          'res' : 0,
    #          'xrng': bb[i]['x_range'],  
    #          'yrng': bb[i]['y_range'],  
    #          'zrng': bb[i]['z_range'],
    #          'id'  : bb[i], 
    #          'mask': A[i]
    #         } for i in range(len(bb)) if con[i]]
    #    with ThreadPool(num_threads) as tp:
    #        out = tp.map(toolbox.getMaskData, di)
    #        sys.stdout.flush() #DEBUG
    #        print(ch) ##DEBUG
    #        sys.stdout.flush() #DEBUG
    #    ChanList.append(np.asarray(out))
    #cubes = np.asarray(ChanList)

    ## For getting bounding box around centroid of annotation
    #for ch in CHAN_NAMES:
    #    di = [{
    #          'rem': rem,
    #          'ch_rsc':
    #            ChannelResource(ch,COLL_NAME,EXP_NAME,'image',datatype='uint8'),
    #          'ch'  : ch,
    #          'res' : 0,
    #          'xrng': [Bcon[i][2] - bf[2], Bcon[i][2] + bf[2] + 1], 
    #          'yrng': [Bcon[i][1] - bf[1], Bcon[i][1] + bf[1] + 1], 
    #          'zrng': [Bcon[i][0] - bf[0], Bcon[i][0] + bf[0] + 1], 
    #         } for i in range(len(Bcon))]
    #    with ThreadPool(num_threads) as tp:
    #        out = tp.map(toolbox.getCube, di)
    #        print(ch) ##DEBUG
    #        sys.stdout.flush() #DEBUG
    #    ChanList.append(np.asarray(out))
    #cubes = np.asarray(ChanList)

    ## for getting all regardles of near boundary
    for ch in CHAN_NAMES:
        di = [{
              'rem': rem,
              'ch_rsc':
                ChannelResource(ch,COLL_NAME,EXP_NAME,'image',datatype='uint8'),
              'ch'  : ch,
              'res' : 0,
              'xrng': [max([Bglobal[i][2] - bf[2], 0]), min([Bglobal[i][2] + bf[2] + 1, ex['x']])], 
              'yrng': [max([Bglobal[i][1] - bf[1], 0]), min([Bglobal[i][1] + bf[1] + 1, ex['y']])], 
              'zrng': [max([Bglobal[i][0] - bf[0], 0]), min([Bglobal[i][0] + bf[0] + 1, ex['z']])] 
             } for i in range(len(Bglobal))]
        with ThreadPool(num_threads) as tp:
            out = tp.map(toolbox.getCube, di)
            print(ch) ##DEBUG
            sys.stdout.flush() #DEBUG
        ChanList.append(np.asarray(out))
    cubes = ChanList
    loc = np.asarray(Bglobal)

    return(cubes, loc)
## END main

def testMain():
    COLL_NAME      = 'collman' 
    EXP_NAME       = 'collman15v2' 
    COORD_FRAME    = 'collman_collman15v2'
    CONFIG_FILE    = 'config.ini'
    OUTPUT         = 'fmaxTest20171214.csv'

    CHAN_NAMES = ['Synapsin647', 'VGluT1_647']
    #CHAN_NAMES = ['DAPI1st', 'DAPI2nd', 'DAPI3rd', 'GABA488', 'GAD647',
    #        'gephyrin594', 'GS594', 'MBP488', 'NR1594', 'PSD95_488',
    #        'Synapsin647', 'VGluT1_647']
    #CHAN_NAMES = ['synapsin', 'PSDr'] 

    cubes, loc  = main(COLL_NAME, EXP_NAME, COORD_FRAME,
         CHAN_NAMES=CHAN_NAMES, num_threads = 6, CONFIG_FILE= 'config.ini')

    Fmaxb = toolbox.Fmaxb(cubes)
    #F0 = toolbox.F0(cubes)
    #Fmax = toolbox.Fmax(cubes)
    #toolbox.mainOUT(Fmax, CHAN_NAMES, OUTPUT)
    #toolbox.toh5(EXP_NAME, OUTPUT + '.h5', CHAN_NAMES, cubes, loc, Fmax)
    return(cubes, loc, Fmaxb)
## End testMain   


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 
        'Get volume normalized F0 values from annotation id regions in the BOSS ')
    parser.add_argument('-C', help='Valid collection id',
            type = str, metavar='C', default='collman')
    parser.add_argument('-E', help='Valid experiment id',
            type = str, metavar='E', default='collman15v2')
    parser.add_argument('-F', help='valid coordinate frame', 
            type = str, metavar='F', default='collman_collman15v2')
    #toolbox.toh5(EXP_NAME, OUTPUT + '.h5', CHAN_NAMES, cubes, loc, F0)
    parser.add_argument('-O', help='output filename',
            type = str, metavar='O', required=True,
            default = 'output')
    parser.add_argument('--con', help='user config file for BOSS'
            'authentication', type = str, metavar='con', required=True)
    
    args = parser.parse_args()

    COLL_NAME      = args.C
    EXP_NAME       = args.E
    COORD_FRAME    = args.F
    OUTPUT         = args.O
    CONFIG_FILE    = args.con

    #rem = BossRemote(CONFIG_FILE)
    #CHAN_NAMES = rem.list_channels(COLL_NAME, EXP_NAME) 

    ##collman15v2 channels
    CHAN_NAMES = ['DAPI1st', 'DAPI2nd', 'DAPI3rd', 'GABA488', 'GAD647',
                  'gephyrin594', 'GS594', 'MBP488', 'NR1594', 'PSD95_488',
                  'Synapsin647', 'VGluT1_647']

    ##collman14v2 channels
    #CHAN_NAMES = ['bIIItubulin', 'DAPI_2nd', 'DAPI_3rd', 
    #              'GABA', 'GAD2', 'VGAT', 'gephyrin', 
    #              'NR1', 'VGluT1',  'synapsin', 'PSDr'] 

    cubes, loc = main(COLL_NAME, EXP_NAME, COORD_FRAME, 
                       CHAN_NAMES=CHAN_NAMES, 
                       num_threads = 6, CONFIG_FILE= CONFIG_FILE)

    F0 = toolbox.F0(cubes) 
    Fmax = toolbox.Fmax(cubes)
    toolbox.mainOUT(Fmax, CHAN_NAMES, OUTPUT)
    idx = np.argsort([3,2,1])
    toolbox.mainOUT(np.transpose(loc[:,idx]), ['x','y','z'], "locations_"+OUTPUT)
    #toolbox.toh5(EXP_NAME, OUTPUT + '.h5', CHAN_NAMES, cubes, loc, F0)
    print('Done!')

