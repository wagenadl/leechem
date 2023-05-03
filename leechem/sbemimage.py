#!/usr/bin/python3

import numpy as np
import urllib.request
import cv2

def _intify(x):
    if type(x)==np.array:
        return x.astype(int)
    else:
        return int(x)
    

class ImageDB:
    '''IMAGEDB - Connection to the SBEM image database online
'''
    def __init__(self, address="http://leechem.caltech.edu:9092"):
        self.address = address
        info = str(self.fetch("info"), "utf-8")
        lines = info.split("\n")
        self.info = {}
        for line in lines:
            kv = line.split("=", 1)
            if len(kv)==2:
                k = kv[0].strip()
                v = kv[1].strip()
                if v.startswith('"'):
                    self.info[k] = v[1:-1]
                elif '.' in v:
                    self.info[k] = self._safefloat(v)
                else:
                    self.info[k] = self._safeint(v)
                    
    def fetch(self, req):
        '''FETCH - Fetch data from server using public API
    FETCH(req), where REQ is a request for the server, returns the server's
    response. For instance, REQ could be “help” or “tile/a2/z2000/y35/x11”.
    The result is BYTES.
'''
        with urllib.request.urlopen(f'{self.address}/{req}') as response:
            return response.read()

    def _safefloat(self, s):
        try:
            return float(s)
        except ValueError:
            return f"<{s}>"

    def _safeint(self, s):
        try:
            return int(s)
        except ValueError:
            return f"<{s}>"

    def tile(self, ix, iy, iz, a=0, b=0):
        '''TILE - Retrieve a 512x512 tile
    Returns a 512x512-pixel tile from the aligned volume at scale
    A (0...8), slice IZ (0...9603), position IX, IY.

    Positions are scale-dependent. At A=0, IX runs from approximately
    2 to 65, and IY runs from approximately 21 to 208, although the
    range is less near the ends of the stack. At A=8, IX = IY = 0 is
    the only tile.

    Optionally, B may be included to average over 2^B slices in Z. In
    that case, IZ runs from 0 to 9603/2^B. If B is nonzero, only
    specific combinations of A and B are supported, corresponding to
    approximately cubic voxels:

        A   B   xy-reso  z-reso (nm)
       --- ---  -------  -----------
        4   1      88      100
        5   2     176      200
        6   3     352      400
        7   4     704      800
        8   5    1408     1600
'''        
        data = self.fetch(f"tile/A{a}/B{b}/Z{iz}/Y{iy}/X{ix}")
        return cv2.imdecode(np.frombuffer(data, np.uint8),
                            cv2.IMREAD_UNCHANGED)

    def roi(self, x, y, iz, w, h, a=0, b=0):
        '''ROI - Retrieve an arbitrarily sized rectangle
    Returns a WxH-pixel image with top-left at (X,Y) from slice IZ.

    Positions are scale-dependent. At A=0, X runs from approximately
    1,000 to 33,000, and Y runs from approximately 10,000 to 106,000, 
    although the range is less near the ends of the stack.

    Optionally, B may be included to average over 2^B slices in IZ. In
    that case, IZ runs from 0 to 9603/2^B. If B is nonzero, only
    specific combinations of A and B are supported, see TILE. 
'''
        data = self.fetch(f"roi_pix/A{a}/B{b}/Z{iz}/Y{y}/X{x}/W{w}/H{h}")
        return cv2.imdecode(np.frombuffer(data, np.uint8),
                            cv2.IMREAD_UNCHANGED)

    def pixtoum(self, x, y, iz, a=0, b=0):
        '''PIXTOUM - Convert pixel position to global µm

        (x, y, z) = PIXTOUM(x, y, iz, a) converts pixel position
        (X,Y) in the slice at IZ to global position in microns,
        suitable for comparison to the tracing database.

        All arguments may be numpy arrays, but if so, they must be
        compatible within each dimension (X with A; Y with A; IZ with B). 

        Restrictions on A and B from ROI do not apply here.
        '''
        sxy = 2**a
        sz = 2**b
        dz = (sz-1)/2/sz # mean displacement of slices in "B" stack
        x = (sxy*x) * self.info['dx']
        y = (sxy*y) * self.info['dy']
        z = (sz*(iz+dz)) * self.info['dz']
        return x, y, z

    def tiletoum(self, dx, dy, ix, iy, iz, a=0, b=0):
        '''TILETOUM - Convert pixel position in tile to global µm
        (x, y, z) = TILETOUM(dx, dy, ix, iy, iz, a)
        converts pixel position (DX,DY) in the tile specified by
        (IX, IY, IZ, A, B) to global position in microns, suitable 
        for comparison to the tracing database. 

        All arguments may be numpy arrays, but if so, they must be
        compatible within each dimension (DX with IX and A; DY with IY
        and A; IZ with B). 

        Restrictions on A and B from TILE do not apply here.
        '''
        return self.pixtoum(dx+512*ix, dy+512*iy, iz, a, b)

    def umtopix(self, x, y, z, a=0, b=0):
        '''UMTOPIX - Convert global µm to pixel position at given scale.
     
        (X, Y, IZ) = UMTOPIX(x, y, z, a, b) converts the given global
        position (X, Y, Z), measured in micrometers as in the tracing
        databse, to pixel positions at the given magnification
        scale. The returned (X, Y, IZ) can be used with ROI,
        assuming (A, B) obey ROI's restrictions.

        All arguments may be numpy arrays, but if so, they must be
        compatible within each dimension (X with A; Y with A; Z with B).

        '''
        sxy = 2**a
        sz = 2**b
        x = (x/self.info['dx']) / sxy
        y = (y/self.info['dy']) / sxy
        z = (z/self.info['dz']) / sz
        return _intify(x), _intify(y), _intify(z)

    def umtotile(self, x, y, z, a=0, b=0):
        '''UMTOTILE - Convert global µm to tile position at given scale.
     
        (DX, DY, IX, IY, IZ) = UMTOTILE(x, y, z, a, b) converts the given
        global position (X, Y, Z), measured in micrometers as in the
        tracing databse, to pixel positions and tile identifiers at the 
        given magnification scale. The returned (IX, IY, IZ) can be used
        with TILE, assuming (A, B) obey TILE's restrictions.

        All arguments may be numpy arrays, but if so, they must be
        compatible within each dimension (X with A; Y with A; Z with B). 
        '''
        x, y, iz = self.umtopix(x, y, z, a, b)
        dx = x % 512
        dy = y % 512
        ix = x // 512
        iy = y // 512
        return dx, dy, ix, iy, iz

    def pixelsize(self, a):
        '''PIXELSIZE - Return size of pixels at given magnification
        PIXELSIZE(a) returns the pixel size (in µm) at the given A-scale.
        '''
        return self.info['dx'] * 2**a

    def slicethickness(self, b):
        '''SLICETHICKNESS - Return thickness of slice at given magnification
        SLICETHICKNESS(b) returns the thickness (in µm) of the
        combined slice at the given B-scale.
        '''
        return self.info['dz'] * 2**b
    
    def voxelsize(self, a, b=0):
        '''VOXELSIZE - Return size of voxels at given magnification
        dx,dy,dz = VOXELSIZE(a,b) returns the 3d-voxel size (in µm)
        at the given A and B-scales.
        '''
        dx = self.info['dx'] * 2**a
        dy = self.info['dy'] * 2**a
        dz = self.info['dz'] * 2**b
        return dx, dy, dz
