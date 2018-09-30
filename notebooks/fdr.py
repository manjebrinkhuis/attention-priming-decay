#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 13:46:11 2016

@author: manje
"""
import numpy as np


def fdr(img, q=.05):
    """ """

    zscores = np.reshape(img, img.size)
    V = len(zscores)
    cV = 1.0    # overwritten
    cV = np.sum([1.0/i for i, z in enumerate(zscores)])
    idx = np.argsort(zscores)
    thrline = idx*q/(V*cV)
    thr = np.max(zscores[zscores <= thrline])

    return thr
