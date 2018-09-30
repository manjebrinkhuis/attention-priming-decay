#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 17:54:26 2017

@author: manje
"""

from nipype.algorithms.modelgen import spm_hrf
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gammaln
from scipy.ndimage import convolve
from roi_selection import Cluster, Sphere
import nibabel as nib
import os

# coordinates from Kristjansson 2007 (roi = (index) # MNI)
clr_lFEF = (62, 66, 62)     # -34, 6, 52
clr_rFEF = (29, 62, 61)     # 32, -2, 50
clr_lIPS = (58, 32, 60)     # -26, -62, 48
clr_rIPS = (25, 39, 65)     # 40, -48, 58
loc_lFEF = (61, 57, 63)     # -32, -12, 54
loc_rFEF = (31, 59, 64)     # 28, -8, 56
loc_lIPS = (60, 33, 56)     # -30, -60, 40
loc_rIPS = (33, 30, 60)     # 24, -66, 48

# Color
color_kristjansson = [
    clr_lFEF, clr_rFEF,
    clr_lIPS, clr_rIPS]

# Location
location_kristjansson = [
    loc_lFEF, loc_rFEF,
    loc_lIPS, loc_rIPS]

# Create ROIs in MNI space
root = "/mnt/data/FMRI/Easter"
base = "loc_lag1"

# Files
mask_fname = os.path.join(root, "masks", "mask.nii.gz")
onset_fname = os.path.join(root, "masks", "mask.nii.gz")
lh_fname = os.path.join(root, "masks", "lh_cortex.nii.gz")
rh_fname = os.path.join(root, "masks", "rh_cortex.nii.gz")
contrast_fname = os.path.join(root, "contrasts", base+".nii.gz")
tmp_fname = "/masked_".join(os.path.split(contrast_fname))
mni_fname = os.path.join(root, "standard", 'MNI152_T1_2mm.nii.gz')

# Mask data
mask = nib.load(mask_fname).get_data()

# Peaks from lag 1 color priming
img = nib.load(contrast_fname)
masked = img.get_data() * mask
img = nib.Nifti1Image(masked, img.affine)
nib.save(img, tmp_fname)

# Cluster
cluster = Cluster(tmp_fname, thresh=-2.3, _min=True, minextent=32)
coords = cluster.get_coords()

# Masks
sph = Sphere(radius=8)
cortex = nib.load(lh_fname).get_data()
cortex += nib.load(rh_fname).get_data()
cortex = convolve(cortex, np.ones((3,3,3)), mode='nearest') > 0
stim_on = nib.load(onset_fname).get_data()

# A lot of masks combined.
super_mask = mask.copy()
super_mask *= cortex.astype(bool)
super_mask *= stim_on.astype(bool)
    
masks = []
for i, (x,y,z) in enumerate(coords):

    # To integers
    x, y, z = int(x), int(y), int(z)
    
    # Create sphere selection in MNI
    x_slice = slice(x-sph.radius, x+sph.radius+1)
    y_slice = slice(y-sph.radius, y+sph.radius+1)
    z_slice = slice(z-sph.radius, z+sph.radius+1)

    # Sphere in MNI
    sphere = np.zeros(mask.shape, dtype=bool)
    sphere[x_slice, y_slice, z_slice] = sph.box.astype(bool)
    sphere *= super_mask.astype(bool)

    # Create nifti image.
    img = nib.Nifti1Image(
            sphere.astype(int),
            nib.load(mni_fname).affine
            )
    
    # Save as nifti.
    roi_fname = os.path.join(root, "rois", "MNI",
                             "%s_cluster_%02d.nii.gz" % (base, i))
    nib.save(img, roi_fname)
    masks.append(roi_fname)

# Convert ROIs from MNI to subject space, per subject.

# Select timeseries within ROI and average across voxels (spatial)

# Concatenate Scans per subject

# Create design matrix:
# 4 regressors for stimulus onset at four locations
# 4 regressors for stimulus repetition; color and location.
# Convolve regressors with canonical HRF.

# Add regressors for scan, accounting for baseline.

