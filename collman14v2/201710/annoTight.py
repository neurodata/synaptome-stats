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
         CHAN_NAMES=None, ANNO_NAME=None, num_threads = 4, CONFIG_FILE= 'config.ini'):

    bf = [5,108,108] # in z,y,x order

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
    anno_res = ChannelResource(ANNO_NAME, COLL_NAME, EXP_NAME, 'annotation', datatype='uint64')
    
    ex = {'x': cfr.x_stop, 'y': cfr.y_stop, 'z': cfr.z_stop}
    
    blocks = block_compute(0,ex['x'],0,ex['y'],0,ex['z'],
               origin = (0,0,0), block_size = (512, 512, 16))

    rid = []
    for b in blocks:
      rid = rid + rem.get_ids_in_region(anno_res, 0, b[0], b[1], b[2], [0,1])

    #u = np.unique(np.asarray(rid[0:4])) ## DEBUG
    u = np.unique(np.asarray(rid))

    ## bounding box for annotation_i
    bb = [rem.get_bounding_box(anno_res, 0,ui, 'tight') for ui in u]

    for i in range(len(bb)):
        bb[i]["id"] = u[i]

    A = [np.uint8(rem.get_cutout(
          anno_res, 0, bb[i]["x_range"], 
          bb[i]["y_range"], bb[i]["z_range"], 
          id_list = [bb[i]['id']]) != bb[i]['id'])
        for i in range(len(bb))] 

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
    conID = []
    con = [np.asarray(mm[j] and MM[j]) for j in range(len(mm))]
    for i in range(len(Bglobal)):
        if con[i]: 
            Bcon.append(Bglobal[i])
            conID.append(bb[i]['id'])

    print(Bcon[0])
    print(conID[0])
    AA = [np.uint8(rem.get_cutout(
           anno_res, 0, 
           [Bcon[i][2] - bf[2], Bcon[i][2] + bf[2]+1],
           [Bcon[i][1] - bf[1], Bcon[i][1] + bf[1]+1],
           [Bcon[i][0] - bf[0], Bcon[i][0] + bf[0]+1],
           id_list = [conID[i]]) == conID[i]) for i in range(len(Bcon))] 
       
    ChanList = []

    ## For getting masked bounding box around centroid of annotation
    for ch in CHAN_NAMES:
        di = [{
              'rem': rem,
              'ch_rsc':
                ChannelResource(ch,COLL_NAME,EXP_NAME,'image',datatype='uint8'),
              'ch'  : ch,
              'res' : 0,
              'xrng': [Bcon[i][2] - bf[2], Bcon[i][2] + bf[2] + 1], 
              'yrng': [Bcon[i][1] - bf[1], Bcon[i][1] + bf[1] + 1], 
              'zrng': [Bcon[i][0] - bf[0], Bcon[i][0] + bf[0] + 1], 
              'mask': AA[i]
             } for i in range(len(Bcon))]
        with ThreadPool(num_threads) as tp:
            out = tp.map(toolbox.getMaskDataT, di)
            print(ch) ## STATUS
            sys.stdout.flush() ## STATUS
        ChanList.append(np.asarray(out))
    cubes = np.asarray(ChanList)
    loc = np.asarray(Bcon)
    nonzeros = [np.sum(ai) for ai in AA]
    return(cubes, loc, nonzeros)
### END main() 


def testMain():
    COLL_NAME      = 'collman' 
    EXP_NAME       = 'collman14v2' 
    ANNO_NAME      = 'annotation' 
    COORD_FRAME    = 'collman_collman14v2'
    CONFIG_FILE    = 'config.ini'
    OUTPUT         = 'testCubes.csv'

    #CHAN_NAMES = ['Synapsin647']
    #CHAN_NAMES = ['DAPI1st', 'DAPI2nd', 'DAPI3rd', 'GABA488', 'GAD647',
    #        'gephyrin594', 'GS594', 'MBP488', 'NR1594', 'PSD95_488',
    #        'Synapsin647', 'VGluT1_647']
    CHAN_NAMES = ['synapsin', 'PSDr'] 

    cubes, loc, nonzeros = main(COLL_NAME, EXP_NAME, COORD_FRAME,
         CHAN_NAMES=CHAN_NAMES, ANNO_NAME = ANNO_NAME, num_threads = 6, CONFIG_FILE= 'config.ini')
    
    F0 = toolbox.F0(cubes)

    #F0w = np.divide(F0, nonzeros)
    toolbox.mainOUT(F0, CHAN_NAMES, OUTPUT)
    toolbox.toh5(EXP_NAME, OUTPUT + '.h5', CHAN_NAMES, cubes, loc, F0, nonzeros)
    
    #import pdb; pdb.set_trace() ## DEBUG
    return(cubes, loc, F0, nonzeros)
## END testMain()


def testPlot(x):
    sns.heatmap(x)
    plt.show()
    return(None)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 
        'Get volume normalized F0 values from annotation id regions in the BOSS ')
    parser.add_argument('-C', help='Valid collection id',
            type = str, metavar='C', default='collman')
    parser.add_argument('-E', help='Valid experiment id',
            type = str, metavar='E', default='collman15v2')
    parser.add_argument('-A', help='Annotation channel name',
            type = str, metavar='A', default='annotation')
    parser.add_argument('-F', help='valid coordinate frame', 
            type = str, metavar='F', default='collman_collman15v2')
    parser.add_argument('-O', help='output filename',
            type = str, metavar='O', required=True,
            default = 'output')
    parser.add_argument('--con', help='user config file for BOSS'
            'authentication', type = str, metavar='con', required=True)
    
    args = parser.parse_args()

    COLL_NAME      = args.C
    EXP_NAME       = args.E
    ANNO_NAME      = args.A
    COORD_FRAME    = args.F
    OUTPUT         = args.O
    CONFIG_FILE    = args.con

    #rem = BossRemote(CONFIG_FILE)
    #CHAN_NAMES = rem.list_channels(COLL_NAME, EXP_NAME) 
    #CHAN_NAMES = ['DAPI1st', 'DAPI2nd', 'DAPI3rd', 'GABA488', 'GAD647',
    #              'gephyrin594', 'GS594', 'MBP488', 'NR1594', 'PSD95_488',
    #              'Synapsin647', 'VGluT1_647']

    ##collman14v2 channels
    CHAN_NAMES = ['bIIItubulin', 'DAPI_2nd', 'DAPI_3rd', 
                  'GABA', 'GAD2', 'VGAT', 'gephyrin', 
                  'NR1', 'VGluT1',  'synapsin', 'PSDr'] 

    cubes, loc, nonzeros = main(
                       COLL_NAME, EXP_NAME, COORD_FRAME, 
                       CHAN_NAMES=CHAN_NAMES, 
                       ANNO_NAME = ANNO_NAME,
                       num_threads = 6, CONFIG_FILE= CONFIG_FILE)

    F0 = toolbox.F0(cubes) 
    toolbox.mainOUT(F0, CHAN_NAMES, OUTPUT)
    toolbox.toh5(EXP_NAME, OUTPUT + '.h5', CHAN_NAMES, cubes, loc, F0)
    print('Done!')

