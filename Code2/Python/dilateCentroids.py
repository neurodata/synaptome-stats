#!/usr/bin/env python3
###
###  
### Jesse Leigh Patsolic For NeuroData
### 2017 <jpatsol1@jhu.edu>
### started from a file made by Alex Baden 
###
### S.D.G 
#
#import grequests # for async requests, conflicts with requests somehow
import itertools
import configparser
import csv
import morton
import sys
import os
import numpy as np 
import requests
import argparse
import time
import tifffile

from numpy import genfromtxt
from numpy import linalg as LA
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import cpu_count
from PIL import Image 

from intern.remote.boss import BossRemote
from intern.resource.boss.resource import * 
from requests import HTTPError 

rmt = BossRemote('neurodata.cfg')

def _3dby8toRGB ( indata ):
    """Convert a numpy array of 3d, 8-bit data to 32bit RGB"""

    rgbdata = np.zeros ( indata.shape[0:2], dtype=np.uint32)

    # rgbdata = np.uint32(0xFF000000) + np.left_shift(np.uint32(indata[:,:,0]),16) + np.left_shift(np.uint32(indata[:,:,1]),8) + np.uint32(indata[:,:,2])
    rgbdata =  np.left_shift(np.uint32(indata[:,:,2]),16) + np.left_shift(np.uint32(indata[:,:,1]),8) + np.uint32(indata[:,:,0])
    return rgbdata

def _RGBto3dby8 ( indata ):
    """Convert a numpy array of 32bit RGB to 3d, 8-bit data"""

    _3ddata = np.zeros ( [indata.shape[0],indata.shape[1],3], dtype=np.uint8)

    _3ddata[:,:,0] = np.uint8(indata&0x000000FF) 
    _3ddata[:,:,1] = np.uint8(np.right_shift(indata&0x0000FF00,8)) 
    _3ddata[:,:,2] = np.uint8(np.right_shift(indata&0x00FF0000,16)) 

    return _3ddata

def read_image(filepath):
    im = Image.open(filepath).convert("RGB")
    imarr = np.asarray(im)
    outslice = _3dby8toRGB(imarr)
    outarr = np.empty((*outslice.shape, 1), dtype=np.uint32)
    outarr[:,:,0] = outslice
    return outarr.astype(np.uint64) 


class SynaptomesBoss:
    def __init__(self, collection, experiment, channel, resolution):
        self.collection = collection 
        self.experiment = experiment
        self.channel = channel
        self.res = resolution

        # Note: An annotation channel must have a source channel. For now we use a blank channel called "image".
        #chan_setup = ChannelResource(self.channel, self.collection, self.experiment, 'annotation', '', sources=['image'], datatype='uint64')
        chan_setup = ChannelResource(self.channel, self.collection, self.experiment, 'image', '', sources=['image'], datatype='uint8')
        try:
            self.chan_actual = rmt.get_project(chan_setup)
        except HTTPError:
            self.chan_actual = rmt.create_project(chan_setup)
        perms = ['read', 'add', 'update', 'add_volumetric_data', 'read_volumetric_data', 'delete_volumetric_data']
        rmt.add_permissions('synaptome', self.chan_actual, perms)
    def post(self, data, x_rng, y_rng, z_rng):
        try:
            rmt.create_cutout(self.chan_actual, self.res, x_rng, y_rng, z_rng, data)
        except Exception as e:
            print("Error: ")
            print(e)


def makeTestImage(filename):
    arr = np.zeros((2048, 2048), dtype=np.uint32)
    # arr = np.random.randint(0, high=1024, size=(2048,2048), dtype=np.uint32)
    val1 = np.random.randint(0, high=1024, dtype=np.uint32)
    arr[500:700, 500:700] = val1 
    val2 = np.random.randint(0, high=1024, dtype=np.uint32)
    arr[1000:1400, 1000:1200] = val2 
    rgbarr = _RGBto3dby8( arr )
    im = Image.fromarray(rgbarr)
    im.save(filename)


def l2ball(BF):
    a = np.zeros((2 * BF + 1, 2 * BF + 1, 2 * BF + 1))
    for i in range(a.shape[0]):
        for j in range(a.shape[1]):
            for k in range(a.shape[2]):
                if LA.norm(np.asarray([k - BF, j - BF, i - BF])) <= BF:
                    a[k,j,i] = 1
                    
    out = np.asarray(a, dtype = np.uint8)
    return(out)

def main():
    parser = argparse.ArgumentParser('Post rectangular annotations stack to the boss')
    parser.add_argument('locations', type=str, help='centroids in x y z csv file')
    parser.add_argument('base_path', type=str, help='Directory where image stacks are to be saved (e.g. "/data/images/"')
    parser.add_argument('base_filename', type=str, help='Base filename with z values appended "fname_z<p:4>"')
    parser.add_argument('--buf', type=int, nargs=3, help='buffer in x y z', default=[5,5,5])
    parser.add_argument('--dim', type=int, nargs=3, help='dimensions in x y z')
    parser.add_argument('--res', type=int, default=0, help='Resolution at which we post data.')
    parser.add_argument('--inftyball', action='store', type=bool, default=False, help='bool to use infinity balls.')
    parser.add_argument('--l2ball', action='store', type=bool, default=False, help='bool to use l2 balls, only uses buf[0].')

    args = parser.parse_args()

    base_path = args.base_path
    base_fname = args.base_filename

    dimx = args.dim[0]
    dimy = args.dim[1]
    dimz = args.dim[2]

    L = genfromtxt(args.locations,  delimiter=',', skip_header = 1, dtype='int32').tolist()
    
    BFx = args.buf[0] 
    BFy = args.buf[1]
    BFz = args.buf[2]

    x_rng = [[x[0] - BFx, x[0] + BFx + 1] for x in L] 
    y_rng = [[y[1] - BFy, y[1] + BFy + 1] for y in L] 
    z_rng = [[z[2] - BFz, z[2] + BFz + 1] for z in L]

    inftyball = np.ones((2 * BFz + 1, 2 * BFy + 1, 2 * BFx + 1), dtype = np.uint8)

    if args.inftyball:
        anno = np.asarray(255*inftyball, dtype = np.uint8)
        ## NUMPY array must be in Z Y X order for the BOSS
    
        z = np.zeros((dimz,dimy,dimx), dtype = np.uint8)
    
        for i in range(len(L)):
            z[z_rng[i][0]:z_rng[i][1],y_rng[i][0]:y_rng[i][1],x_rng[i][0]:x_rng[i][1]] = anno

        fout = os.path.join(base_path, base_fname)

        for j in range(z.shape[0]):
            tmp = fout + "_z{}.tif".format(str(j).zfill(len(str(z.shape[0]))))
            tifffile.imsave(tmp, np.asarray(z[j,:,:], dtype = np.uint8))

        tifffile.imsave(fout + "_stack.tif", np.asarray(z, dtype = np.uint8))
        print("Test image is here: {}".format(base_path))
        sys.exit(0)

    if args.l2ball:
        anno = 255*l2ball(BFx) 
        z = np.zeros((dimz,dimy,dimx), dtype = np.uint8)
    
        for i in range(len(L)):
            rep = anno == 255
            z[z_rng[i][0]:z_rng[i][1],y_rng[i][0]:y_rng[i][1],x_rng[i][0]:x_rng[i][1]][rep] = 255
            
        fout = os.path.join(base_path, base_fname)

        for j in range(z.shape[0]):
            tmp = fout + "_z{}.tif".format(str(j).zfill(len(str(z.shape[0]))))
            tifffile.imsave(tmp, np.asarray(z[j,:,:], dtype = np.uint8))

        tifffile.imsave(fout + "_stack.tif", np.asarray(z, dtype = np.uint8))
        print("Test image is here: {}".format(base_path))
        sys.exit(0)

    
    return(None)


if __name__ == '__main__':
    main()




## python3 dilateCentroids.py ~/tmps/uptmp/Ex12R75_buff5.csv ~/tmps/uptmp/ Synapsin1_2_synapses --buf 5 5 5 --dim 5491 4749 35 --inftyball True
