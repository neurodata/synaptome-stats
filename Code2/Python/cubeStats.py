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




def cube(D):

    CONFIG_FILE = D['config']

    BFx = D['bf'] ## FIX THIS
    BFy = D['bf'] ## FIX THIS
    BFz = D['bf'] ## FIX THIS
    x_rng = [D['x']-BFx, D['x']+BFx+1] 
    y_rng = [D['y']-BFy, D['y']+BFy+1] 
    z_rng = [D['z']-BFz, D['z']+BFz+1]
    COLL_NAME = D['COLL_NAME']
    EXP_NAME = D['EXP_NAME']
    CHAN_NAMES = D['CHAN_NAMES']
    res = D['res']

    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    TOKEN = config['Default']['token']
    boss_url = ''.join( ( config['Default']['protocol'],'://',config['Default']['host'],'/v1/' ) )
    #print(boss_url)
    #'https://api.boss.neurodata.io/v1/'
    
    #intern
    rem2 = BossRemote(CONFIG_FILE)

    cubeCh = []

    for ch in CHAN_NAMES:
        ch_rsc = ChannelResource(ch,COLL_NAME,EXP_NAME,'image',datatype='uint16')
        ### Get cube and transpose to X, Y, Z order
        data = rem2.get_cutout(ch_rsc, res, x_rng, y_rng, z_rng)
        
        #data = np.transpose(data, axes=(2,1,0)
        #outArray = np.reshape(data, 
        #        (z_rng[1] - z_rng[0], 
        #         y_rng[1] - y_rng[0], 
        #         x_rng[1] - x_rng[0]), 
        #        order = 'C')

        cubeCh.append(data)

    out = np.asarray(cubeCh)
    print([D['x'], D['y'], D['z']]) ##DEBUG
    sys.stdout.flush() ##DEBUG
    return(out)


def test():
    COLL_NAME  = 'weiler14'
    EXP_NAME  = 'Ex2R18C1'
    COORD_FRAME  = 'weiler14_Ex2R18C1' 
    LOCATIONS_FILE = 'locationTest.csv'
    BF = 5
    CONFIG_FILE = 'config.ini'

    CHAN_NAMES= ['Synapsin1_3', 'Synaptopodin_6', 'vGluT1_3'],
    LOCATIONS = [[619, 242, 21], [1150, 1995, 21]]

    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    TOKEN = config['Default']['token']
    boss_url = ''.join( ( config['Default']['protocol'],'://',config['Default']['host'],'/v1/' ) )
    #print(boss_url)
    #'https://api.boss.neurodata.io/v1/'
    
    #intern
    rem1 = BossRemote(CONFIG_FILE)
    #ch_list = rem1.list_channels(COLL_NAME, EXP_NAME)
    #out = cube(D)

    cubes = main(COLL_NAME, EXP_NAME, COORD_FRAME, LOCATIONS, BF = BF,
            num_threads = 4, CONFIG_FILE=CONFIG_FILE)

    #out = np.transpose(cubes, axes=(4,3,2,1,0))

    #print("Cube data has been transposed to x y z channel sample") ## DEBUG
    print('Saving data')
    mainOUT(COORD_FRAME, "jtest_output", cubes, LOCATIONS)
    print("Done")
    return(cubes)



def main(COLL_NAME, EXP_NAME, COORD_FRAME, LOCATIONS, BF = 5,
         CHAN_NAMES=None, num_threads = 4, CONFIG_FILE= 'config.ini'):

    L = LOCATIONS
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    TOKEN = config['Default']['token']
    boss_url = ''.join( ( config['Default']['protocol'],'://',config['Default']['host'],'/v1/' ) )
    #print(boss_url)
    #'https://api.boss.neurodata.io/v1/'
    
    #intern
    rem1 = BossRemote(CONFIG_FILE)

    res=0

    coll_rsc = CollectionResource(COLL_NAME)
    exp_rsc = ExperimentResource(COLL_NAME,EXP_NAME,COORD_FRAME)

    
    if CHAN_NAMES is None:
        #CHAN_NAMES = rem1.list_channels(COLL_NAME, EXP_NAME)
        CHAN_NAMES = ['DAPI_1', 'DAPI_2', 'DAPI_3', 'DAPI_4', 'DAPI_5a', 'DAPI_5b',
                      'DAPI_6', 'DAPI_7', 'GABAARa1_7', 'GAD2_4', 'Gephyrin_1',
                      'GFP_5b', 'GluR1_5a', 'GluR2_6', 'GluR4_7', 'NR2A_2',
                      'NR2B_4', 'PSD25_2', 'PV25_1', 'Synapsin1_3',
                      'Synaptopodin_6', 'vGAT_3', 'vGluT1_3', 'vGluT2_2', 'YFP_1']

    di = [] 

    for a in range(len(L)):
        di.append({'config': CONFIG_FILE, 
            'COLL_NAME': COLL_NAME, 
            'EXP_NAME': EXP_NAME,
            'COORD_FRAME': COORD_FRAME,
            'CHAN_NAMES': CHAN_NAMES,
            'res': 0, 
            'x': int(L[a][0]), 
            'y': int(L[a][1]), 
            'z': int(L[a][2]),
            'bf': int(BF)
            })

    with ThreadPool(num_threads) as tp:
        out = tp.map(cube, di)

    outb = np.asarray(out)
    return(outb)

def testH5():
    a = np.arange(100).reshape(2,5,10)
    f = h5py.File("hdf5TEST_2_5_10.h5", 'w')
    f.create_dataset("test", data = a)
    f.close()

    
def mainOUT(TOKEN,OUTPUT, out, L):
    h5f0OUT = h5py.File(OUTPUT+".h5", 'w')
    h5f0OUT.create_dataset(TOKEN + "_cubes", data = out)
    h5f0OUT.create_dataset("Locations", data = np.transpose(L))
    h5f0OUT.close()
    return(None)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 
        'Download synapse cubes using the BOSS')
    parser.add_argument('-C', help='Valid collection id',
            type = str, metavar='C', default='weiler14')
    parser.add_argument('-E', help='Valid experiment id',
            type = str, metavar='E', default='Ex2R18C1')
    parser.add_argument('-F', help='valid coordinate frame', 
            type = str, metavar='F', default='weiler14_Ex2R18C1')
    parser.add_argument('-B', help='integer buffer around center',
            type = str, metavar='B', default=5)
    parser.add_argument('-L', help='csv file of locations '
            'in xyz order', type = str, metavar='L', required=True)
    parser.add_argument('-O', help='output filename',
            type = str, metavar='O', required=True)
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

    LOCATIONS = genfromtxt(LOCATIONS_FILE, delimiter=',', skip_header = 0).tolist()

    cubes = main(COLL_NAME, EXP_NAME, COORD_FRAME, LOCATIONS, BF = BF,
                 CHAN_NAMES=None, num_threads = 12, CONFIG_FILE= CONFIG_FILE)

    print("Cube data is formatted in `x y z channel sample` order.")
    print('Saving data')
    mainOUT(COORD_FRAME, OUTPUT, cubes, LOCATIONS)
    print("Done")


