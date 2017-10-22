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


def getData(di):
    data = di['rem'].get_cutout(di['ch_rsc'], di['res'], di['xrng'],
          di['yrng'] ,di['zrng'])

    out = data[di['mask']]

    return(out)

def getCentroid(box):
    m = np.asarray(box == True) 
    avg = np.int(np.round(np.mean(m,1)))
    return(avg)



def main(COLL_NAME, EXP_NAME, COORD_FRAME,
         CHAN_NAMES=None, num_threads = 4, CONFIG_FILE= 'config.ini'):

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
    anno_res = ChannelResource('annotation', 'collman', 'collman15v2', 'annotation', datatype='uint64')
    
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
            
    Cen = np.asarray([getCentroid(bi) for bi in bb])

    #A = [rem.get_cutout(anno_res, 0, bb[i]["x_range"], bb[i]["y_range"],
    #                    bb[i]["z_range"], id_list = [bb[i]['id']]) != 0 for i
    #                    in range(len(bb))]

    #if CHAN_NAMES is None:
    #            CHAN_NAMES = ['DAPI1st', 'DAPI2nd', 'DAPI3rd',
    #    	      'GABA488', 'GAD647', 'gephyrin594', 'GS594', 'MBP488',
    #                  'NR1594', 'PSD95_488', 'Synapsin647', 'VGluT1_647']

    #ChanList = []
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
    #         } for i in range(len(bb))]
   
    #    with ThreadPool(num_threads) as tp:
    #        out = tp.map(getData, di)
    #        print(ch) ##DEBUG
    #        sys.stdout.flush() #DEBUG

    #    ChanList.append(np.asarray(out))

    #anno = np.asarray(ChanList)

    #F0 = np.asarray([[sum(a)/len(a) for a in anno[i]] for i in range(anno.shape[0])])
    #return(F0)
    return(Cen)

def testMain():
    COLL_NAME      = 'collman' 
    EXP_NAME       = 'collman15v2' 
    COORD_FRAME    = 'collman_collman15v2'
    CONFIG_FILE    = 'config.ini'

    CHAN_NAMES = ['PSD95_488', 'Synapsin647', 'VGluT1_647']

    anno = main(COLL_NAME, EXP_NAME, COORD_FRAME,
         CHAN_NAMES=CHAN_NAMES, num_threads = 4, CONFIG_FILE= 'config.ini')

    F0 = np.asarray([[sum(a)/len(a) for a in anno[i]] for i in range(anno.shape[0])])

    return(F0)
   

def mainOUT(F0, head, outfile):

    with open(outfile, 'w') as f1:
        wt = csv.writer(f1)
        wt.writerow(head)
        wt.writerows(np.transpose(F0))

    return(None)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 
        'Get volume normalized F0 values from annotation id regions in the BOSS ')
    parser.add_argument('-C', help='Valid collection id',
            type = str, metavar='C', default='collman')
    parser.add_argument('-E', help='Valid experiment id',
            type = str, metavar='E', default='collman15v2')
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
    COORD_FRAME    = args.F
    OUTPUT         = args.O
    CONFIG_FILE    = args.con

    rem = BossRemote(CONFIG_FILE)

    #CHAN_NAMES = rem.list_channels(COLL_NAME, EXP_NAME) 
    CHAN_NAMES = ['DAPI1st', 'DAPI2nd', 'DAPI3rd', 'GABA488', 'GAD647',
                  'gephyrin594', 'GS594', 'MBP488', 'NR1594', 'PSD95_488',
                  'Synapsin647', 'VGluT1_647']

    F0 = main(COLL_NAME, EXP_NAME, COORD_FRAME, 
                       CHAN_NAMES=CHAN_NAMES, 
                       num_threads = 4, CONFIG_FILE= CONFIG_FILE)

    #mainOUT(F0, CHAN_NAMES, OUTPUT)
    mainOUT(F0, ['z', 'y', 'x'], OUTPUT)
    print('Done!')

