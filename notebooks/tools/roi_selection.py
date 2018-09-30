#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 17:21:42 2017

@author: manje
"""

import numpy as np
import subprocess as sp
from nibabel import Nifti1Image


def dict2opts(
    dct, 
    doubledash=True, 
    change_underscores=True, 
    assign_sign='='):
    """ dictionary to options """
    
    opts = []
    for kw in dct.keys():
        # use long or short options?
        dash = '--' if doubledash else '-'
        # remove first underscore if it is present
        opt = kw[1:] if kw[0] == '_' else kw
        if change_underscores:
            # change remaining underscores to single dash
            opt = opt.replace('_', '-')
        if not type(dct[kw]) == bool:
            opts.append('{0}{1}{2}{3}'.format(dash, opt, 
                                              assign_sign, 
                                              dct[kw]))
        elif dct[kw] is True:
            opts.append('{0}{1}'.format(dash, opt))

    return opts


def get_radii(coords):
    """ Radii of x,y,z arrays in array. Distance from (0,0,0). """
    return [(x**2+y**2+z**2)**.5 for x,y,z in coords]


class Sphere():
    """"""
    def __init__(self, radius=8):
        
        self.radius = int(radius)
        self.diameter = radius*2+1
        self.box = np.ones((radius*2+1,)*3)
        self.rad_box = np.ones((radius*2+1,)*3)
        
        coords = [(int(x), int(y), int(z)) 
                  for x in range(self.diameter) 
                  for y in range(self.diameter) 
                  for z in range(self.diameter)]

        coords = np.array(coords, dtype=int)
        distances = get_radii(coords - self.radius)
        
        for i, (x,y,z) in enumerate(coords):
            if distances[i] > (self.radius):
                self.box[x,y,z] = 0
            else:
                self.box[x,y,z] = 1

        self.rad_box[x,y,z] = distances[i]


class Cluster():
    """ Use FSL to create clusters. """

    def __init__(self, in_file, out_file='index.txt', 
                 thresh=1.96, **kwargs):
        self.in_file = in_file
        self.out_file = out_file
        self.thresh = thresh
        self.cmd = ['cluster', '--in='+in_file, '--thresh='+str(thresh)]
        self.cmd += dict2opts(kwargs)
        self.data = None

    def run(self):
        pipe = sp.Popen(' '.join(self.cmd), shell=True, stdout=sp.PIPE)
        lines = pipe.stdout.read().split('\n')
        hdr = lines[0].split('\t')
        data = [tuple(line.split('\t')) for 
                line in lines[1:-1]]
        self.data = np.array(data, 
                             dtype=[(name, 'float') for name in hdr])
    
    def get_coords(self):
        
        if self.data is None:
            self.run()
            
        coords = self.data[["MAX X (vox)",
                            "MAX Y (vox)",
                            "MAX Z (vox)"]]
        
        return coords
    
            
        
        
        