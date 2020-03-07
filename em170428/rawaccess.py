#!/usr/bin/python3
# -*- coding: utf-8 -*-

import cv2
import sys
import diced
import numpy as np
import subprocess

dcd = diced.DicedStore("dvid://http://127.0.0.1", port=8000)
REPONAME = 'raw170428'
try:
    repo = dcd.open_repo(REPONAME)
except diced.DicedException as exc:
    print(exc)
    print('Constructing repo')
    aa = ['dvid', 'repos', 'new', REPONAME, 'Raw data for the 170428 SBEM run']
    res = subprocess.run(aa, check=True)
    repo = dcd.open_repo(REPONAME)

REPO = repo.uuid
IMAGEWIDTH = 17100
QIMAGEWIDTH = IMAGEWIDTH//4
QQIMAGEWIDTH = QIMAGEWIDTH//4

def rawFilename(r, m, s):
    '''RAWFILENAME - Path to original raw file
    RAWFILENAME(r, m, s) returns the pathname of the original raw image for
    the given run, montage, slice.'''
    IFNBASE="/lsi1/push/170428-SBEM/Run%i/Montage_%03i/Run%i_OnPoint_%04i.tif"
    return IFNBASE % (r, m, r, s)

def stackName(r, m):
    '''STACKNAME - DVID name of image stack for a montage.
    STACKNAME(r, m) returns the DVID name of the image stack for the given
    run, montage'''
    return "r%im%i" % (r, m)

def loadImage(ifn):
    '''LOADIMAGE - Load raw image
    LOADIMAGE(ifn) loads the named raw image.'''
    return cv2.imread(ifn, cv2.IMREAD_ANYDEPTH + cv2.IMREAD_GRAYSCALE)

def stretchImage(img, stretch=None):
    '''STRETCHIMAGE - Convert 16 to 8 bit and stretch contrast
    STRETCHIMAGE(img) returns an 8-bit version of the given image.
    STRETCHIMAGE(img, stretch) also stretches the contrast, dropping 
    STRETCH percent of dark and light pixels. In this case, output image
    uses pixel values 1 through 255 so that 0 can serve as "transparent."'''
    if stretch:
        N = img.size
        ilo = int(.01*stretch*N)
        ihi = int((1-.01*stretch)*N)
        hst = cv2.calcHist([(img//256).astype('uint8')], 
                           [0], None, [256], [0,256])
        cum = np.cumsum(hst)
        vlo = 256*np.argmax(cum>=.01*stretch*cum[-1])
        vhi = 256*np.argmax(cum>=(1-.01*stretch)*cum[-1]) - 1
        vlo -= 256
        if vlo<0:
            vlo = 0
        if vhi==0:
            vhi = 65535
            vhi += 256
        if vhi > 65535:
            vhi = 65535
        nrm = np.array([255./(vhi-vlo)]).astype('float32')
        img -= vlo
        img *= 255. / (vhi - vlo)
        img[img<1] = 1
        img[img>255] = 255
        img = img.astype('uint8')
    else:
        img = img.astype('uint8')
    return img

def blurImage(img, rad=3):
    '''BLURIMAGE - Gaussian blur
    BLURIMAGE(img, radius) applies a Gaussian blur with given radius
    and default sigma to the image.'''
    return cv2.GaussianBlur(img, (rad,rad), 0)

def scaleImage(img, fac=4):
    '''SCALEIMAGE - Downscale an image
    SCALEIMAGE(img, fac) down scales an image by the given integer factor 
    (default: 4) using full area averaging.'''
    S = img.shape
    Ky = S[0]//fac
    Kx = S[1]//fac
    if fac*Ky<S[0] or fac*Kx<S[1]:
        img = img[0:fac*Ky,0:fac*Kx]
    return cv2.resize(img, None, None, 1/fac, 1/fac, cv2.INTER_AREA)

def pushImage(img, irun, imont, islic):
    '''PUSHIMAGE - Store image in DVID server
    PUSHIMAGE(img, irun, imont, islic) pushes the given image as the
    raw image for the given run, montage, slice. Also pushes 1/4² and 1/16²
    resolution versions.'''
    stk = stackName(irun,imont)
    arr = repo.get_array(stk)
    S = img.shape
    arr[islic,0:S[0],0:S[1]] = np.reshape(img, [1, S[0], S[1]])
    ldd = repo.get_array('r%iloaded' % irun)
    ldd[0,islic,imont] = np.array([[[1]]], dtype='uint8')

    img = scaleImage(img)
    stk += 'q'
    arr = repo.get_array(stk)
    S = img.shape
    arr[islic,0:S[0],0:S[1]] = np.reshape(img, [1, S[0], S[1]])
    ldd = repo.get_array('r%iqloaded' % irun)
    ldd[0,islic,imont] = np.array([[[1]]], dtype='uint8')

    img = scaleImage(img)
    stk += 'q'
    arr = repo.get_array(stk)
    S = img.shape
    arr[islic,0:S[0],0:S[1]] = np.reshape(img, [1, S[0], S[1]])
    ldd = repo.get_array('r%iqqloaded' % irun)
    ldd[0,islic,imont] = np.array([[[1]]], dtype='uint8')
    

def contained(irun, imont, islic):
    '''CONTAINED - Test if image already uploaded
    CONTAINED(irun, imont, islic) returns True if and only if an image
    has already been uploaded for the given run, montage, slice.'''
    try:
        ldd = repo.get_array('r%iloaded' % irun)
        return ldd[0,islic,imont] > 0
    except:
        return False
    
def dvidRepoCommand(args):
    '''DVIDREPOCOMMAND - Run a DVID command directly
    DVIDREPOCOMMAND(args) runs a "dvid repo" command directly on our repo.'''
    aa = ['dvid', 'repo', REPO]
    aa += args
    res = subprocess.run(aa, stderr=subprocess.PIPE,
                         stdout=subprocess.PIPE, check=True)
    return (res.stdout, res.stderr)

def getRawImage(irun, imont, islic, rect=None):
    '''GETRAWIMAGE - Download raw image from DVID
    GETRAWIMAGE(irun, imont, islic) downloads the entire raw image for the
    given run, montage, slice.
    GETRAWIMAGE(irun, imont, islic, (x0, y0, w, h)) specifies that only
    the given rectangle be downloaded.
    Result is a numpy array.'''
    stk = stackName(irun,imont) + 'q'
    arr = repo.get_array(stk)
    if rect is None:
      return arr[islic, 0:QIMAGEWIDTH, 0:QIMAGEWIDTH]
    else:
      return arr[islic, rect[1]:rect[1]+rect[3], rect[0]:rect[0]+rect[2]]
  
def getQuarter(irun, imont, islic):
    '''GETQUARTER - Download raw image at 1/4² resolution
    GETQUARTER(irun, imont, islic) downloads the raw image for the given
    run, montage, slice at 1/4² resolution.'''
    stk = stackName(irun,imont) + 'q'
    arr = repo.get_array(stk)
    return arr[islic, 0:QIMAGEWIDTH, 0:QIMAGEWIDTH]

def getQQuarter(irun, imont, islic):
    '''GETQQUARTER - Download raw image at 1/16² resolution
    GETQQUARTER(irun, imont, islic) downloads the raw image for the given
    run, montage, slice at 1/16² resolution.'''
    stk = stackName(irun,imont) + 'qq'
    arr = repo.get_array(stk)
    return arr[islic, 0:QQIMAGEWIDTH, 0:QQIMAGEWIDTH]
  
if __name__ == '__main__':
  # This example pushes a single raw image to the server
  # Call as "rawaccess R M S".
  irun = int(sys.argv[1])
  imont = int(sys.argv[2])
  islic = int(sys.argv[3])
  
  img = loadImage(rawFilename(irun, imont, islic))
  if img is None:
    raise ValueError('Image not found')
  
  img = blurImage(img, 3)
  img = stretchImage(img, .2)
  pushImage(img, irun, imont, islic)
  
