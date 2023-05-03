#!/usr/bin/python3

import numpy as np
import h5py as h5
import math
import sys
import os.path
from . import webaccess

TILESIZE = 512
IMAGESIZE = TILESIZE * 33
IMAGERAD = IMAGESIZE / 2
RAWSIZE = 17100
MARGIN = (RAWSIZE - IMAGESIZE)//2

class BackTransform:
    def __init__(self, x0,y0,bx,cx,by,cy):
        '''BACKTRANSFORM - Class to contain back-style perspective transforms
        BACKTRANSFORM(x0, y0, bx, cx, by, cy) prepares a transform object
        that can be used to convert model coordinates to source image coords.'''
        self.x0 = x0
        self.y0 = y0
        self.bx = bx
        self.cx = cx
        self.by = by
        self.cy = cy
    def modelToSource(self, x, y):
        '''MODELTOSOURCE - Map model coordinates to source image coordinates
        (x, y) = MODELTOSOURCE(x, y) transforms the point (x,y) from model
        coordinates (that is, the space where Maria has been working) to
        source image coordinates (pixels in the 17100 x 17100 space of the 
        raw images.)'''
        xr = x - self.x0
        yr = y - self.y0
        xi = (xr-IMAGERAD)/IMAGERAD
        eta = (yr-IMAGERAD)/IMAGERAD
        return (xr - self.bx*eta - self.cx*xi*eta + MARGIN,
                yr - self.by*xi - self.cy*xi*eta + MARGIN)
    def sourceContains(self, x, y):
        (x1,y1) = self.modelToSource(x, y)
        return x1>=0 and x1<RAWSIZE and y1>=0 and y1<RAWSIZE

class RunInfo:
    f = None
    nruns = None
    z2r = None
    zz0 = None
    
    def __init__(self, ifn=None):
        '''RUNINFO - Class to wrap back position data for the 170428 run
        pp = RUNINFO(ifn) creates a new instance of the RUNINFO
        class, loading the data from the given hdf5 file.
        pp = RUNINFO() loads the data from the 'positionsummary.h5' 
        file in the em170428 directory.'''
        ifn = webaccess.ensurefile(ifn, "positionsummary.h5")
        self.f = h5.File(ifn, 'r')

    BLEND = 1024
    
    def blender(self, x):
        print('blend', x, self.BLEND)
        if x<=0:
            return 0
        elif x>=self.BLEND:
            return 1
        else:
            return .5 - .5*math.cos(x*math.pi/self.BLEND)

    def runCount(self):
        '''RUNCOUNT - Number of runs in dataset
        n = pp.RUNCOUNT() returns the number of runs in the dataset.'''
        if self.nruns is None:
            self.nruns = 0
            for k in self.f['runs'].keys():
                r = int(k[1:])
                if r>self.nruns:
                    self.nruns = r
        return self.nruns

    def _ensure_z2r(self):
        if self.z2r is None:
            nr = self.runCount()
            ns = 0
            for r0 in range(nr):
                ns += self.sliceCount(r0+1)
            self.z2r = np.zeros(ns, dtype='int32')
            self.zz0 = np.zeros(nr, dtype='int32')
            z0 = 0
            for r0 in range(nr):
                s = self.sliceCount(r0+1)
                self.z2r[z0:z0+s] = r0 + 1
                self.zz0[r0] = z0
                z0 += s

    def zToRunSlice(self, z):
        '''ZTORUNSLICE - Convert Z coordinate to run and slice
        (r,s) = ZTORUNSLICE(z) returns the run and slice number for a given
        Z position. Remember that runs count from one but slices from zero.
        If Z is a list or a numpy array, the resulting R and S are of the
        same type.'''
        if type(z)==list:
            r = []
            s = []
            for z1 in z:
                (r1,s1) = self.zToRunSlice(z1)
                r.append(r1)
                s.append(s1)
            return (r,s)
        elif type(z)==np.ndarray:
            N = z.size
            zz = z.reshape((N,))
            r = 0*zz
            s = 0*zz
            for k in range(N):
                (r[k], s[k]) = self.zToRunSlice(zz[k])
            r = r.reshape(z.shape)
            s = s.reshape(z.shape)
            return (r,s)

        self._ensure_z2r()
        r = self.z2r[z]
        s = z - self.zz0[r-1]
        return (r,s)

    def runSliceToZ(self, r, s):
        self._ensure_z2r()
        z0 = self.zz0[r-1]
        return z0 + s
    
    def sliceCount(self, irun):
        '''SLICECOUNT - Number of slices in a run
        n = pp.SLICECOUNT(r) returns the number of slices in run R.
        Remember that runs are counted from 1.'''
        r = self.f['runs']['R%i' % irun]
        return r['x0'].shape[0]

    def montageCount(self, irun):
        '''MONTAGECOUNT - Number of montages in a run
        n = pp.MONTAGECOUNT(r) returns the number of montages in run R.
        Remember that runs are counted from 1.'''
        r = self.f['runs']['R%i' % irun]
        return r['x0'].shape[1]

    def columnCount(self, irun):
        '''COLUMNCOUNT - Number of columns in a run
        n = pp.COLUMNCOUNT(r) returns the number of columns in run R.
        Remember that runs are counted from 1.'''
        m = self.montageCount(irun)
        if m<=3:
            return 1
        else:
            return 2

    def rowCount(self, irun):
        '''ROWCOUNT - Number of rows in a run
        n = pp.ROWCOUNT(r) returns the number of rows in run R.
        Remember that runs are counted from 1.'''
        m = self.montageCount(irun)
        if m<=3:
            return m
        else:
            return m//2

    def transform(self, irun, imontage, islice):
        '''TRANSFORM - Return transform for a given source image
        xf = TRANSFORM(irun, imontage, islice) returns the transformation
        for the given source image. 
        XF is of type BackTransform and can be used to transform points in
        model space to source image space.'''
        r = self.f['runs']['R%i' % irun]
        x0 = r['x0'][islice,imontage]
        y0 = r['y0'][islice,imontage]
        bx = r['bx'][islice,imontage]
        cx = r['cx'][islice,imontage]
        by = r['by'][islice,imontage]
        cy = r['cy'][islice,imontage]
        # See E&R p.1510 on why I cannot return an affine or perspective
        return BackTransform(x0, y0, bx, cx, by, cy)

    def tileCenter(self, irun, imontage, islice):
        '''TILECENTER - Find center of a bigtile
           (x,y,z) = TILECENTER(r, m, s) returns the pixel location of the
           bigtile specified by run, montage, and slice.
           (A "bigtile" is 33x33 small 512x512 tiles.)'''
        r = self.f['runs']['R%i' % irun]
        x = r['x0'][islice,imontage] + IMAGERAD
        y = r['y0'][islice,imontage] + IMAGERAD
        self._ensure_z2r()
        z = self.zz0[irun-1] + islice
        return (x,y,z)
        
    def findXYZ(self, xyz):
        '''FINDXYZ - Find montages that contain a certain point
        (r, s, mm) = FINDXYZ(xyz) where XYZ is an (x,y,z) triplet in model
        space returns the run R, slice S, and details on montages that
        identify the raw images that went into building the pixel at (x,y,z).
        Montage info is returned as a dict of montage numbers
        to a (x,y,alpha) tuple, where (X,Y) are pixel coordinates in the
        raw source image and ALPHA is the blend used. (ALPHA = 1 if only one
        image affected the final pixel, otherwise the ALPHAs of all returned
        montages add up to one.)'''
        (r, s) = self.zToRunSlice(xyz[2])
        M = self.montageCount(r)
        mm = {}
        xy = []
        salpha = 0
        for m in range(M):
            xf = self.transform(r, m, s)
            (x,y) = xf.modelToSource(xyz[0], xyz[1])
            if x>=0 and x<RAWSIZE and y>=0 and y<RAWSIZE:
                mm[m] = (x, y, 1)
        if len(mm)>1:
            alph = {}
            for m in mm.keys(): 
                xf = self.transform(r, m, s)
                (x,y) = xf.modelToSource(xyz[0], xyz[1]) 
                x1 = x - MARGIN
                y1 = y - MARGIN
                al = 1
                if x1 < self.BLEND:
                    al *= self.blender(x1)
                elif x1 > IMAGESIZE - self.BLEND:
                    al *= self.blender(IMAGESIZE - x1)
                if y1 < self.BLEND:
                    al *= self.blender(y1)
                elif y1 > IMAGESIZE - self.BLEND:
                    al *= self.blender(IMAGESIZE - y1)
                alph[m] = al
                salpha += al
            for m in mm.keys():
                xya = mm[m]
                mm[m] = (xya[0], xya[1], alph[m] / (salpha+1e-99))
        return (r, s, mm)
    
    def rawToTile(self, xy):
        '''RAWTOTILE - Find tile coordinates for raw image coordinates
        (XY, xy) = RAWTOTILE(xy) find tile numbers and pixel coordinates
        in smoothed unaligned jpegs given raw coordinates.
        Recall that the raw images measured 17100x17100 pixels, that a margin
        of 102 pixels was chopped off on all four sides, leaving a grid
        of 33x33 tiles each measuring 512x512 pixels.'''
        x0 = xy[0] - MARGIN
        y0 = xy[1] - MARGIN
        XY = (int(x0)//TILESIZE,
              int(y0)//TILESIZE)
        xy = (x0 - TILESIZE*XY[0],
              y0 - TILESIZE*XY[1])
        return (XY, xy)

if __name__=='__main__':
    import matplotlib.pyplot as plt
    import cv2
    pp = BackPositions('positionsummary.h5')
    (r,s,mm) = pp.findXYZ((36088,116776,3657))
    print(r,s,mm)
    m = []
    for m1 in mm.keys():
        m = m1
    (XY, xy) = pp.rawToTile(mm[m])
    print(XY, xy)
    
    img = cv2.imread('/md0/dw/170428-SBEM/unaligned/R%03i/M%03i/S%04i/X%02iY%02i.tif' % (r, m, s, XY[0], XY[1]))

    plt.imshow(img)
    plt.plot(xy[0], xy[1],'r*')

    ######
    plt.figure()
    (r,s,mm) = pp.findXYZ((38408,96904,1477))
    print(r,s,mm)
    m = []
    for m1 in mm.keys():
        m = m1
    (XY, xy) = pp.rawToTile(mm[m])
    print(XY, xy)
    img = cv2.imread('/md0/dw/170428-SBEM/unaligned/R%03i/M%03i/S%04i/X%02iY%02i.tif' % (r, m, s, XY[0], XY[1]))

    plt.imshow(img)
    plt.plot(xy[0], xy[1],'r*')
    
