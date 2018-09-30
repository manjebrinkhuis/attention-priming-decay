#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 13:54:06 2016

@author: manje
"""

import nibabel as nib
from clusters import FSLCluster
import os
from nipype.interfaces import fsl
import numpy as np

neg = False
thresh = 3.1

localizer = os.path.join('localizer_L3_fixedfx_reg2',
                         '{loc}', '_fixedflameo0',
                         'zstat1.nii.gz')

target = os.path.join('decay_final2_fe',
                      '{loc}', '_fixedflameo0',
                      'zstat1.nii.gz')

intersected = 'intersected_{loc}.nii.gz'
stats = 'stats_{loc}.txt'
out_index = 'out_index_{loc}.nii.gz'

# Intersect target data and localizer
for con in [0,1,2,3]:

    continue

    print con

    loc = nib.load(localizer.format(loc=con))
    tgt = nib.load(target.format(loc=con))

    if neg:
        loc *= -1
        tgt *= -1

    mask = (loc.get_data() > thresh)
    masked = tgt.get_data() * mask.astype(float)
    masked2 = masked > thresh

    img = nib.Nifti1Image(masked, tgt.get_affine())
    nib.save(img, intersected.format(loc=con))

    cluster = FSLCluster(in_file=intersected.format(loc=con),
                         thresh=thresh,
                         out_file=stats.format(loc=con),
                         minextent=10,
                         oindex=out_index.format(loc=con))
    cluster.run()

    masked2

# %% Convert ROIs to subject space
reg = 'registration/{0}'
struct = 'T1_hires.nii.gz'
inplane = 'inplane_brain.nii.gz'
field = 'field.nii.gz'
inv_field = 'field_inverse.nii.gz'
premat = 'premat.mat'
subjects = ['001',
            '002',
            '003',
            '004',
            '005',
            '006']

fname = 'out_index_{loc}.nii.gz'
ROI_locs = [fname.format(loc=i) for i in range(4)]

print 'Converting ROIs to subject space...'
for sub in subjects:
    print '- subject', sub
    path = reg.format(sub)
    inv_warp = fsl.InvWarp(warp=os.path.join(path, field),
                           inverse_warp=os.path.join(path, inv_field),
                           reference=os.path.join(path, struct),
                           output_type='NIFTI_GZ')
    #inv_warp.run()

    postmat = np.linalg.inv(np.loadtxt(os.path.join(path, 'premat.mat')))
    np.savetxt(os.path.join(path, 'postmat.mat'), postmat, fmt='%.8f')
    out_path = os.path.join('rois', sub)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    for roi in ROI_locs:
        print '\t- ROI', roi
        warp = fsl.ApplyWarp(field_file=os.path.join(path, inv_field),
                             ref_file=os.path.join(path, inplane),
                             in_file=os.path.join('rois', roi),
                             out_file=os.path.join(out_path, roi),
                             postmat=os.path.join(path, 'postmat.mat'),
                             interp='nn')

        #warp.run()

