#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 12:12:37 2017

@author: manje
"""

import pandas as pd
import os
import numpy as np
from priming import(priming_per_condition, priming_per_location,
                    priming, decay, repeat_array, array2dummies,
                    create_priming_model)
import matplotlib.pyplot as plt
from hrf import create_regressor


# Files
root = "/mnt/data/FMRI/Easter"
combined_fname = os.path.join(root, "combined.csv")
data = pd.read_csv(combined_fname)
subs = sorted(list(set(data['sub'])))

for s in subs:
    print "=== sub %d ===" % s
    
    sub = data[data["sub"] == s]
    idx = sub.first_valid_index()
    block = sub["block"]
    block = block - block[idx]
    
    M = sub[["loc1", "loc2", "loc3", "loc4"] +
            ["loc_lag1", "loc_lag2",
             "clr_lag1", "clr_lag2"]]
        
    blocks = sorted(list(set(block)))
    for b in blocks:
        M["reg%d" % b] = np.array(block == b).astype(int)

    M = np.array(M).T
    
    clusters1 = ["clr_kristjansson_cluster_00_normalized",
                 "clr_kristjansson_cluster_01_normalized",
                 "clr_kristjansson_cluster_02_normalized",
                 "clr_kristjansson_cluster_03_normalized"]
    names1 = ["lFEF", "rFEF", "lIPS", "rIPS"]    
    Y1 = sub[clusters1]
    
    clusters2 = ["loc_kristjansson_cluster_00_normalized",
                 "loc_kristjansson_cluster_01_normalized",
                 "loc_kristjansson_cluster_02_normalized",
                 "loc_kristjansson_cluster_03_normalized"]
    Y2 = sub[clusters2]
    names2 = ["lFEF", "rFEF", "lIPS", "rIPS"]
    
    clusters3 = ["clr_lag1_cluster_00_normalized",
                 "clr_lag1_cluster_01_normalized",
                 "clr_lag1_cluster_02_normalized",
                 "clr_lag1_cluster_03_normalized",
                 "clr_lag1_cluster_04_normalized",
                 "clr_lag1_cluster_05_normalized"]
    Y3 = sub[clusters3]
    names3 = ["rSPL", "lSPL", "rFEF", "lFEF", "rLOC", "lLOC"]
        
    Y4 = sub[["loc_lag1_cluster_00_normalized",
              "loc_lag1_cluster_01_normalized",
              "loc_lag1_cluster_02_normalized",
              "loc_lag1_cluster_03_normalized",
              "loc_lag1_cluster_04_normalized",
              "loc_lag1_cluster_05_normalized",
              "loc_lag1_cluster_06_normalized"]]
    names4 = ["rSPL", "lSPL", "lFEF", "lFEF/SMA", "rFEF", "rLOC", "rV1"]
    
    Ys = [Y1, Y2, Y3, Y4]
    names = [names1, names2, names3, names4]
    
    for Yi, Y in enumerate(Ys):
        plt.figure(Yi)
        bs = glm(np.array(Y).T, M).T        
        contrasts = [
#                     [1, -.33, -.33, -.33] + [0] * (len(M)-4),
#                     [-.33, 1, -.33, -.33] + [0] * (len(M)-4),
#                     [-.33, -.33, 1, -.33] + [0] * (len(M)-4),
#                     [-.33, -.33, -.33, 1] + [0] * (len(M)-4),
                     [.25, .25, .25, .25] + [0] * (len(M)-4),
                     [0, 0, 0, 0, 1, 0, 0, 0] + [0] * (len(M)-8),
                     [0, 0, 0, 0, 0, 1, 0, 0] + [0] * (len(M)-8),
                     [0, 0, 0, 0, 0, 0, 1, 0] + [0] * (len(M)-8),
                     [0, 0, 0, 0, 0, 0, 0, 1] + [0] * (len(M)-8)]
    
        for i, y in enumerate(np.array(Y).T):
            zscores = [contrast(y, bs[i], M, contrasts=con)
                       for con in contrasts]
              
            plt.subplot(len(Y.T)/2 + len(Y.T)%2, 2, i+1)
            plt.title(names[Yi][i])
            plt.bar(np.arange(len(zscores)) + 1.0/(len(subs)+2)*s,
                    zscores,
                    width=1.0/(len(subs)+1))
            plt.ylim(-5, 5)

