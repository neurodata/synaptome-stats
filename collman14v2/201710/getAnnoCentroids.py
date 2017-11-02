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
         CHAN_NAME=None, CONFIG_FILE= 'config.ini'):

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
    anno_res = ChannelResource(CHAN_NAME, COLL_NAME, EXP_NAME, 'annotation', datatype='uint64') 
    
    ex = {'x': cfr.x_stop, 'y': cfr.y_stop, 'z': cfr.z_stop}
    
    blocks = block_compute(0,ex['x'],0,ex['y'],0,ex['z'],
               origin = (0,0,0), block_size = (512, 512, 16))

    rid = []
    for b in blocks:
      rid = rid + rem.get_ids_in_region(anno_res, 0, b[0], b[1], b[2], [0,1])

    u = np.unique(np.asarray(rid))

    bb = [rem.get_bounding_box(anno_res, 0,ui, 'tight') for ui in u]

    for i in range(len(bb)):
        bb[i]["id"] = u[i]

    A = [rem.get_cutout(anno_res, 0, bb[i]["x_range"], bb[i]["y_range"],
                        bb[i]["z_range"], id_list = [bb[i]['id']]) != 0 for i
                        in range(len(bb))]

    Bmeans = [np.int32(np.round(np.mean(np.asarray(np.where(A[i] == True)),1))) for i in range(len(A))]
    
    Bglobal = []
    for i in range(len(bb)):
        ad1 = np.asarray([bb[i]['z_range'][0], bb[i]['y_range'][0], bb[i]['x_range'][0]])
        Bglobal.append(Bmeans[i] + ad1)
        
    out = np.asarray(Bglobal) ## in zyx order
    xyz = out[:, [2,1,0]] ## in xyz order

    return(xyz)

def testMain():
    COLL_NAME      = 'collman' 
    EXP_NAME       = 'collman14v2' 
    COORD_FRAME    = 'collman_collman14v2'
    CHAN_NAME      = 'annotation'
    CONFIG_FILE    = 'config.ini'
    OUTPUT         = 'collman14v2_annotationCentroids.csv'

    loc  = main(COLL_NAME, EXP_NAME, COORD_FRAME,
         CHAN_NAME=CHAN_NAME, num_threads = 6, CONFIG_FILE= 'config.ini')

    return(loc)
   


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 
        'Get volume normalized F0 values from annotation id regions in the BOSS ')
    parser.add_argument('-C', help='Valid collection id',
            type = str, metavar='C', default='collman')
    parser.add_argument('-E', help='Valid experiment id',
            type = str, metavar='E', default='collman15v2')
    parser.add_argument('-F', help='valid coordinate frame', 
            type = str, metavar='F', default='collman_collman15v2')
    parser.add_argument('-N', help='valid channel name', 
            type = str, metavar='N', default='annotation')
    parser.add_argument('-O', help='output filename',
            type = str, metavar='O', required=True,
            default = 'output')
    parser.add_argument('--con', help='user config file for BOSS'
            'authentication', type = str, metavar='con', required=True)
    
    args = parser.parse_args()

    COLL_NAME      = args.C
    EXP_NAME       = args.E
    COORD_FRAME    = args.F
    CHAN_NAME      = args.N
    OUTPUT         = args.O
    CONFIG_FILE    = args.con


    loc = main(COLL_NAME, EXP_NAME, COORD_FRAME, CHAN_NAME=CHAN_NAME, 
               CONFIG_FILE= CONFIG_FILE)

    np.savetxt(OUTPUT, loc, delimiter = ',', fmt = "%d", 
               header = "x,y,z")

    print('Done!')
