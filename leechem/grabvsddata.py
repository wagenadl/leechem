#!/usr/bin/python3

import pynwb
import hdf5storage
import numpy as np

src = '/lab/tomina/170428/'

lbLdir = src + 'Localbend/P2L_anatomically/'
lbLtri = [9, 10, 11]
lbRdir = src + 'Localbend/P2R_anatomically/'
lbRtri = [12]
swimdir = src + 'Swim/'
swimtri = [6, 8]
crawldir = src + 'Crawling/'
crawltri = [15, 17]

trialtypes = ['LB-L', 'LB-R', 'Swim', 'Crawl']
dirs = [lbLdir, lbRdir, swimdir, crawldir]
trials = [lbLtri, lbRtri, swimtri, crawltri]

def roiids(mat):
    rois = []
    for a in mat[0]:
        rois.append(a[0])
    return rois

DTR = pynwb.hdmf.common.table.DynamicTableRegion

for typ in range(len(trialtypes)):
    for k in range(len(trials[typ])):
        typlbl = trialtypes[typ]
        dr = dirs[typ]
        tri = trials[typ][k]
        mc = hdf5storage.loadmat(f'/tmp/v2m-{tri}.mat')
        
        ids = DTR(name='roiid',
                  data=roiids(mc['id_vsd']),
                  description='ROI IDs')
        tt = mc['t_vsd'].flatten()
        dFF = np.transpose(mc['dFF_vsd'])
        vsd = pynwb.ophys.RoiResponseSeries(name='VSD',
                                            data=dFF,
                                            timestamps=tt,
                                            rois=ids,
                                            unit='raw')
        
