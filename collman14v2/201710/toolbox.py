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
import h5py
import sys
import os
import itertools
from scipy.stats import multivariate_normal
from functools import partial
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
import csv
import datetime


def getAnnoData(di):
    data = di['rem'].get_cutout(di['ch_rsc'], di['res'], di['xrng'],
          di['yrng'] ,di['zrng'])
    out = np.multiply(data, di['mask'])
    return(out)

def getMaskData(di):
    data = di['rem'].get_cutout(di['ch_rsc'], di['res'], di['xrng'],
          di['yrng'] ,di['zrng'])
    out = np.multiply(data, di['mask'])
    return(out)

def getMaskDataT(di):
    data = di['rem'].get_cutout(di['ch_rsc'], di['res'], di['xrng'],
          di['yrng'] ,di['zrng'])
    out = np.multiply(data, di['mask']).astype(data.dtype)
    return(out)

def getCube(di):
    data = di['rem'].get_cutout(di['ch_rsc'], di['res'], di['xrng'],
          di['yrng'] ,di['zrng'])
    return(data)

def getWCube(di):
    data = di['rem'].get_cutout(di['ch_rsc'], di['res'], di['xrng'],
          di['yrng'] ,di['zrng'])
    y = np.multiply(data, di['w']) 
    return(y)

def getCentroid(box):
    m = np.asarray(box == True) 
    avg = np.int(np.round(np.mean(m,1)))
    return(avg)

def weightCubes(cubes, w):
    c = np.float32(cubes)
    for i in range(c.shape[0]):
        for j in range(c.shape[1]):
            tmp = np.multiply(c[i,j,::], w)
            c[i,j,::] = tmp
    return(c)

def distMat2(bf):
    A = np.reshape(np.array([np.sqrt((i-(bf+1))**2 + (j-(bf+1))**2)
    for i in range(1, 2*bf+2)
    for j in range(1, 2*bf+2)]), (2*bf+1, 2*bf+1))
    A[np.where(A == 0)] = 1
    return(A)

def distMat3(bf):
    A = np.reshape(
        np.array([np.sqrt((i-(bf+1))**2 + (j-(bf+1))**2 + (k-(bf+1))**2)\
        for i in range(1, 2*bf+2)\
        for j in range(1, 2*bf+2)\
        for k in range(1, 2*bf+2)]),(2*bf+1, 2*bf+1, 2*bf+1))
    #A[np.where(A == 0)] = 1
    return(A)

def F0(cubes):
    print(cubes.shape)
    f0 = [[np.sum(cubes[i,j,:,:,:]) for j in range(cubes.shape[1])] for i in range(cubes.shape[0])]
    #f0 = np.sum(cube)
    return(f0)

def oob():

    return(None)

def gaussian3dMask(mean, sd, x,y,z,bf):
    def mvn(x, m = [0,0,0], sigma = [1,1,1]):
            y = multivariate_normal.pdf(x, mean = m, cov = sigma)
            return(y)

    mask = np.array([0 for i in range(x*y*z)], dtype = np.float32).reshape(x,y,z)
    for i in range(x):
      for j in range(y):
        for k in range(z):
          mask[i,j,k] = mvn(x = [i - bf[0], j - bf[1], k - bf[2]])
    return(mask)

def toh5(EXP_NAME, outfile, CHAN_NAMES, cubes, loc, F0 = None, w = None):
    hf = h5py.File(outfile, 'w')
    hf.create_dataset(EXP_NAME + "_cubes", data = cubes)
    hf.create_dataset("Locations", data = np.transpose(loc))
    hf.create_dataset("Channels", data = np.string_(CHAN_NAMES))
    if F0 is not None:
        hf.create_dataset("F0", data = F0)
    if w is not None:
        hf.create_dataset("w", data = w)

    hf.close()
    return(None)

def mainOUT(F0, head, outfile):

    with open(outfile, 'w') as f1:
        wt = csv.writer(f1)
        wt.writerow(head)
        wt.writerows(np.transpose(np.asarray(F0)))

    return(None)




