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
import svgwrite
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

    BF = D['bf']
    x_rng = [D['x']-BF, D['x']+BF+1] 
    y_rng = [D['y']-BF, D['y']+BF+1] 
    z_rng = [D['z']-BF, D['z']+BF+1]
    CHAN_NAME = D['CHAN_NAME']
    COLL_NAME = D['COLL_NAME']
    EXP_NAME = D['EXP_NAME']
    res = D['res']

    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    TOKEN = config['Default']['token']
    boss_url = ''.join( ( config['Default']['protocol'],'://',config['Default']['host'],'/v1/' ) )
    #print(boss_url)
    #'https://api.boss.neurodata.io/v1/'
    
    #intern
    rem2 = BossRemote(CONFIG_FILE)

    ch_rsc = ChannelResource(CHAN_NAME,COLL_NAME,EXP_NAME,'image',datatype='uint16')
    ### Get cube and transpose to X, Y, Z order
    data = np.transpose(rem2.get_cutout(ch_rsc, res, x_rng, y_rng, z_rng), axes=(2,1,0))

    return(data)


def test():
    D = {'config': 'config.ini', 
         'x': 605,
         'y': 603,
         'z': 15,
         'bf': 3,
         'res': 0,
         'CHAN_NAME': 'DAPI_1',
         'COLL_NAME': 'weiler14',
         'EXP_NAME' : 'Ex2R18C1',
         'COORD_FRAME': 'weiler14_Ex2R18C1'}

    out = cube(D)
    print(out)
    return(out)


def main(COLL_NAME, EXP_NAME, COORD_FRAME, LOCATIONS_FILE, BF = 5, 
        num_threads = 4, CONFIG_FILE= 'config.ini'):

    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)
    TOKEN = config['Default']['token']
    boss_url = ''.join( ( config['Default']['protocol'],'://',config['Default']['host'],'/v1/' ) )
    #print(boss_url)
    #'https://api.boss.neurodata.io/v1/'
    
    #intern
    rem1 = BossRemote(CONFIG_FILE)

    L = genfromtxt(LOCATIONS_FILE, delimiter=',', skip_header = 0).tolist()
    res=0

    coll_rsc = CollectionResource(COLL_NAME)
    exp_rsc = ExperimentResource(COLL_NAME,EXP_NAME,COORD_FRAME)

    
    ch_list = rem1.list_channels(COLL_NAME, EXP_NAME)

    L = genfromtxt(LOCATIONS_FILE, delimiter=',', skip_header = 0).tolist()

    di = [] 

    for ch in ch_list:
        for a in range(len(L)):
            di.append({'config': CONFIG_FILE, 
                'COLL_NAME': COLL_NAME, 
                'EXP_NAME': EXP_NAME,
                'COORD_FRAME': COORD_FRAME,
                'CHAN_NAME': ch,
                'res': 0, 
                'x': int(L[a][0]), 
                'y': int(L[a][1]), 
                'z': int(L[a][2]),
                'bf': int(BF)
                })

    with ThreadPool(num_threads) as tp:
        out = tp.map(cube, di)

    return(out)


    
def mainOUT(TOKEN,OUTPUT, out, L):
    h5f0OUT = h5py.File(OUTPUT+".h5", 'w')
    h5f0OUT.create_dataset(TOKEN + "_cubes", data = out)
    h5f0OUT.create_dataset("Locations", data = np.transpose(L))
    h5f0OUT.close()
    return(None)


if __name__ == '__main__':

    COLL_NAME  = 'weiler14'
    EXP_NAME  = 'Ex2R18C1'
    COORD_FRAME  = 'weiler14_Ex2R18C1' 
    LOCATIONS_FILE = 'locationTest.csv'
    BF = 5
    CONFIG_FILE = 'config.ini'

    #print(test())
    cubes = main(COLL_NAME, EXP_NAME, COORD_FRAME, LOCATIONS_FILE, BF = 5,
            num_threads = 4, CONFIG_FILE= 'config.ini')

    L = genfromtxt(LOCATIONS_FILE, delimiter=',', skip_header = 0).tolist()
    mainOUT('weiler14_Ex2R18C1', 'tmp', cubes, L)
    print(cubes[0])









