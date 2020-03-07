#!/usr/bin/python3

import numpy as np
import h5py
import errno
import os

class UCT:
    def __init__(self, h5fn=None):
        '''UCT - Python access to the uCT soma position file
        uct = UCT(h5fn) opens the given "uCT-somata.h5" file.
        uct = UCT() opens the default file in the em170428 directory.
        '''
        if h5fn is None:
            here = os.path.dirname(__file__)
            h5fn = here + '/../data/uCT-somata.h5'
            with h5py.File(h5fn, 'r') as h5:
                self.somapos = h5['somapos']['value'][:] * 2
                self.exitpoint = h5['dpp']['value'][:]
                self.lut = h5['dc']['value'][:]
                self.tails = {}
                tt = h5['tails']['value']
                for k in tt.keys():
                    if k.startswith('_'):
                        k1= int(k[1:])
                        t = tt[k]
                        t = t['value']
                        self.tails[k1]= t[:]
