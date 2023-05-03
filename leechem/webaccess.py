#!/usr/bin/python3

import urllib.request
import shutil
import os
from contextlib import contextmanager
import codecs

@contextmanager
def opentextfile(fn, dflt):
    if fn is not None:
        yield open(fn, mode)
    here = os.path.dirname(__file__)
    localfn = here + "/../data/" + dflt
    if os.path.isfile(localfn):
        yield open(localfn, mode)
    url = "https://leechem.caltech.edu/170428/pydata/" + dflt
    try:
        req = urllib.request.urlopen(url)
        yield codecs.iterdecode(req, "utf-8")
    except urllib.request.HTTPError:
        raise RuntimeError(f"Could not access {url}")
    
        
def ensurefile(fn, dflt):
    if fn is not None:
        return fn
    here = os.path.dirname(__file__)
    localfn = here + "/../data/" + dflt
    if os.path.isfile(localfn):
        return localfn
    
    url = "https://leechem.caltech.edu/170428/pydata/" + dflt
    localfn = "./" + dflt
    try:
        with urllib.request.urlopen(url) as response:
            with open(localfn, 'wb') as fd:
                shutil.copyfileobj(response, fd)
    except urllib.request.HTTPError:
        raise RuntimeError(f"Could not download {url}")
    return localfn
