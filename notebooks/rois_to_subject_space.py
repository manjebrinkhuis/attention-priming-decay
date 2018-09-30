#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 22:34:41 2017

@author: manje
"""

import os
import nibabel as nib
from nipype.interfaces import fsl
import numpy as np

# Convert ROIs from MNI to subject space, per subject.

# Files
root = "/mnt/data/FMRI/Easter"
files = [os.path.join(root, "rois", "MNI", fname) for fname 
         in os.listdir(os.path.join(root, "rois", "MNI"))
         if ".nii.gz" in fname]
mni_fname = os.path.join(root, "standard", 'MNI152_T1_2mm.nii.gz')
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

# Convert to subject space
for sub in subjects:
    print sub
    path = os.path.join(root, "registration", sub)    
    inv_warp = fsl.InvWarp(warp=os.path.join(path, field),
                           inverse_warp=os.path.join(path, inv_field),
                           reference=mni_fname,
                           output_type='NIFTI_GZ')
    inv_warp.run()
    
    postmat = np.linalg.inv(np.loadtxt(os.path.join(path, 'premat.mat')))
    np.savetxt(os.path.join(path, 'postmat.mat'), postmat, fmt='%.8f')
    
    out_path = os.path.join(root, "rois", sub)
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    for f in files:
        _, fname = os.path.split(f)
        print fname
        warp = fsl.ApplyWarp(field_file=os.path.join(path, inv_field),
                             ref_file=os.path.join(path, inplane),
                             in_file=f,
                             out_file=os.path.join(out_path, fname),
                             postmat=os.path.join(path, 'postmat.mat'),
                             interp='nn')

        warp.run()
