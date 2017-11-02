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
import configparser
import requests
import numpy as np
from numpy import genfromtxt
import h5py
import sys
import os
import csv
import morton

def main(locations):

    L = genfromtxt(locations,  delimiter=',', skip_header = 1, dtype='int32').tolist()
    m = morton.Morton(dimensions=3, bits=64)

    Lzo = sorted([m.pack(L[i][0], L[i][1], L[i][2]) for i in range(len(L))])

    loc = [m.unpack(l) for l in Lzo]

    return(loc)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description = 'Sort a list of 3-dimensional points into Morton Z-order')
    parser.add_argument('--locations', type=str, help='Valid path to a csv file with columns in x,y,z order.  Expects a header.')
    parser.add_argument('--output', type=str, help='output csv file name')
    
    args = parser.parse_args()

    locations = args.locations 
    outf = args.output

    out = np.asarray(main(locations), dtype = np.uint32)

    head = 'x,y,z'
    np.savetxt(outf, out, header = head, delimiter=',', fmt = "%d")

    print('Done!')

