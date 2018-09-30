# -*- coding: utf-8 -*-
"""
Created on Wed Jul 27 14:34:49 2016

@author: manje
"""

import nibabel as nib
import os
import numpy as np

# info
roi_dir = 'rois/{sub}'
ROI_names = ['out_index']
subjects = ['001',
            '002',
            '003',
            '004',
            '005',
            '006']

# copes
cope = 'cope{0}.nii.gz'
copes_stim = [cope.format(i) for i in range(1, 5)]
loc_copes_lag1 = [cope.format(i) for i in range(5, 9)]
loc_copes_lag2 = [cope.format(i) for i in range(9, 13)]
clr_copes_lag1 = [cope.format(i) for i in range(13, 17)]
clr_copes_lag2 = [cope.format(i) for i in range(17, 21)]
copes = [copes_stim, loc_copes_lag1, loc_copes_lag2,
         clr_copes_lag1, clr_copes_lag2]

# pes
pe = 'pe{0}.nii.gz'
pes_stim = [pe.format(i) for i in range(1, 9, 2)]
loc_pes_lag1 = [pe.format(i) for i in range(9, 17, 2)]
loc_pes_lag2 = [pe.format(i) for i in range(17, 25, 2)]
clr_pes_lag1 = [pe.format(i) for i in range(25, 33, 2)]
clr_pes_lag2 = [pe.format(i) for i in range(33, 41, 2)]
pes = [pes_stim, loc_pes_lag1, loc_pes_lag2, clr_pes_lag1, clr_pes_lag2]

# pes
cope = 'cope{0}.nii.gz'
loc_pes_stim = [cope.format(i) for i in range(1, 5)]
loc_pes_priming = [cope.format(i) for i in range(5, 9)]
item_copes = [loc_pes_stim, loc_pes_priming]

model = 'decay_items_final_filmgls'
#model = 'decay_loc_final2_filmgls'

img_files = pes

# get means in ROIs
all_data = []
for sub in subjects:
    print 'subject', sub

    # Put scans from session
    # folders in single list
    count_scan = 0
    for ses in ['001', '002']:

        sub_dir = 'ses_{ses}_sub_{sub}'.format(ses=ses, sub=sub)
        path = os.path.join(model, sub_dir)
        scans = [os.path.join(path, d) for d in os.listdir(path)]

        # Go through scans
        for iscan, scan in enumerate(scans):
            print '\tscan', scan

            # For each location,
            for loc in range(4):
                print '\t\tlocation', loc

                # get data
                copes_img = [nib.load(os.path.join(scan, c[loc])) for c in img_files]

                # Other locs.
                other_locs = range(4)
                other_locs.remove(loc)

                # For each ROI
                for roi in ROI_names:
                    print '\t\t\tROI', roi

                    # target mask
                    fname = '{roi}_{loc}.nii.gz'
                    roi_path = os.path.join(roi_dir.format(sub=sub), fname)
                    img = nib.load(roi_path.format(loc=loc, roi=roi))
                    tgt_mask = img.get_data().astype(bool)

                    # distractor mask
                    dis_masks = []
                    for oloc in other_locs:
                        img = nib.load(roi_path.format(loc=oloc, roi=roi))
                        dis_masks.append(img.get_data().astype(bool))

                    tgt_activity = []
                    dis_activity = []
                    for i, img in enumerate(copes_img):
                        data = img.get_data()
                        tgt_activity.append(np.percentile(data[tgt_mask], 90))
                        print 'tgt activity', tgt_activity
                        dis = [np.percentile(data[dm], 90) for dm in dis_masks]
                        print 'dis activity', np.mean(dis), dis
                        dis_activity.append(np.mean(dis))

                    line_tgt = [sub, int(ses), iscan,
                                count_scan, roi, 1] + tgt_activity
                    line_dis = [sub, int(ses), iscan,
                                count_scan, roi, 0] + dis_activity

                    all_data.append(line_tgt)
                    all_data.append(line_dis)

            count_scan += 1

output = [['subject', 'ses', 'scan_ses', 'scan',
           'roi', 'target',
           'stim_on',
           'loc_lag1', 'loc_lag2',
           'clr_lag1', 'clr_lag2']]
output += all_data

fname = 'mixed_model_input_item_copes.csv'
open(fname, 'w').close()
for line in output:
    with open(fname, 'a') as f:
        f.write('\t'.join([str(val) for val in line]) + '\n')
