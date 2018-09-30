#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 16:33:40 2016

@author: manje
"""

import nibabel as nib
from scipy.stats import norm
import matplotlib.pyplot as plt
from matplotlib.colors import colorConverter
import numpy as np
import matplotlib as mpl
from matplotlib import cm
import os

from mpl_toolkits.mplot3d import axes3d


def fdr(img, q=.05):
    """ """

    zscores = np.sort(np.reshape(img, img.size))
    pval = norm.sf(zscores)[::-1]
    V = len(pval)
    idx = np.arange(1, V+1)

    print '='*20

    fig = plt.figure(2)
    ax = fig.add_subplot(111)

    cV = 1.0    # overwritten
    cV = np.sum([1.0/(i+1) for i, p in enumerate(pval)])

    thrline = idx*q/(V*cV)
    select = pval <= thrline
    print 'num voxels:', np.sum(select), 'out of', len(img)

    if np.sum(select) == 0:
        thr = -10
    else:
        thr = np.min(zscores[select])

    ax.plot(idx, pval, '-')
    ax.plot(idx, thrline, '-')
    ax.plot([0, len(zscores)], [0.005, 0.005], '-')

    print 'threshold:', thr

    return thr


def plot_color(ax, color):

    ax.spines['bottom'].set_color(color)
    ax.spines['top'].set_color(color)
    ax.spines['right'].set_color(color)
    ax.spines['left'].set_color(color)

    ax.tick_params(axis='x', colors=color)
    ax.tick_params(axis='y', colors=color)

    ax.yaxis.label.set_color(color)
    ax.xaxis.label.set_color(color)

    ax.title.set_color(color)


class Con2Img():

    def __init__(self,
                 filename,
                 contrasts=[''],
                 fdr_corrected=False,
                 mask_fname=None,
                 background='',
                 ntiles=4,
                 center=60,
                 axis=0,
                 skip=4,
                 thresh=2.3,
                 neg=False,
                 intersect=False,
                 cmaps=['hot'],
                 cmap_back='gray'):

        self.filename = filename
        self.contrasts = contrasts
        self.background = background
        self.cmaps = cmaps
        self.cmap_back = cmap_back
        self.neg = neg
        self.ntiles = ntiles
        self.center = center
        self.axis = axis
        self.skip = skip
        self.thresh = []

        # Contrast image
        self.img_contrasts = []
        for contrast in contrasts:

            img_con = nib.load(contrast).get_data()

            if fdr_corrected:

                include = np.ones(img_con.shape, dtype=bool)
                include *= img_con != 0

                if mask_fname is not None:
                    include *= nib.load(mask_fname).get_data().astype(bool)

                data = img_con[include] * -1 if neg else img_con[include]
                thresh = fdr(data) # * -1 if neg else fdr(data)

            self.thresh.append(thresh)
            mask = img_con > -thresh if neg else img_con < thresh
            img_con[mask] = 0

            self.img_contrasts.append(img_con.astype(float))

        self.intersection = np.zeros(self.img_contrasts[0].shape, dtype=bool)
        if intersect:
            self.intersection = np.invert(self.intersection)
            for con in self.img_contrasts:
                print 'num voxels intersection:', np.sum(self.intersection)
                self.intersection *= con < 0 if neg else con > 0

        # Background image
        self.img_back = nib.load(background).get_data()

    def get_slice(self, index=60, axis=2):

        if axis == 0:
            select = [index, slice(None), slice(None, None, -1)]
            back = self.img_back[select].T
            isect = self.intersection[select].T
            fronts = [img[select].T for img in self.img_contrasts]

            x, y = back.shape
            square = np.zeros((y, y))
            margin = int(abs(y - x))
            nback = square.copy()
            nisect = square.copy()
            nback[margin:, :] = back
            nisect[margin:, :] = isect
            nfronts = [square.copy() for i in range(len(fronts))]

            for i, f in enumerate(fronts):
                nfronts[i][margin:, :] = f

        elif axis == 1:
            select = [slice(None), index, slice(None, None, -1)]
            nback = self.img_back[select].T
            nisect = self.intersection[select].T
            nfronts = [img[select].T for img in self.img_contrasts]

        elif axis == 2:
            select = [slice(None), slice(None, None, -1), index]
            back = self.img_back[select].T
            isect = self.intersection[select].T
            fronts = [img[select].T for img in self.img_contrasts]

            x, y = back.shape
            square = np.zeros((x, x))
            margin = int(abs(y - x)/2)
            nback = square.copy()
            nisect = square.copy()
            nback[:, margin:-margin] = back
            nisect[:, margin:-margin] = isect

            nfronts = [square.copy() for i in range(len(fronts))]
            for i, f in enumerate(fronts):
                nfronts[i][:, margin:-margin] = f

        fronts = nfronts
        back = nback
        isect = nisect

        return back, fronts, isect

    def save_image(self):

        width = self.ntiles * (self.skip+1)
        start = int(self.center - width/2)

        cols = int(self.ntiles / (self.ntiles ** .5))
        rows = int(np.ceil(self.ntiles / float(cols)))

        # Color map
        cmaps = []
        for cmap_name in self.cmaps:
            cmap = plt.imshow(np.ones((16, 16)), cmap=cmap_name).cmap
            cmap._init()
            plt.close()
            if self.neg:
                cmap._lut[::-1] = cmap._lut
            cmap._lut[-4:, -1] = 0
            cmaps.append(cmap)

        indices = range(start, start+width, self.skip+1)

        # Create tiled image
        fig = plt.figure(figsize=((cols)*3 + 3*cols*.3, rows*3),
                         facecolor="#000000")

        for row in range(rows):
            for col in range(cols):

                index = row*(cols)+col
                back, fronts, isect = self.get_slice(indices[index], self.axis)
                plt.subplot(rows, cols, row*(cols)+(col+1))
                back_ax = plt.imshow(back,
                                     cmap=self.cmap_back,
                                     interpolation='nearest')
                back_ax.set_clim(0, 30000)

                axes = []
                for i, front in enumerate(fronts):
                    ax = plt.imshow(front, cmap=cmaps[i],
                                    interpolation='nearest')

                    ax.set_clim(*sorted([self.thresh[i], self.thresh[i]*2]))
                    axes.append(ax)

                # print np.any(isect), cmaps[-1].name
                # isect_ax = plt.imshow(isect.astype(float), cmap=cmaps[-1])

                plt.axis('off')
                plt.title('xyz'[self.axis]+'={0}'.format(indices[index]),
                          color='#ffffff', position=(.8, .8))

        padding = .05
        spacing = 0
        bar_width = .025
        bar_height = 1.0 / rows - padding
        width = .7
        bar_pad = .1
        bottom = .05

        plt.subplots_adjust(left=padding, right=width-padding,
                            top=1-padding, bottom=bottom,
                            hspace=spacing, wspace=spacing)

        plt.axis('off')

        labels = ['z-score lag-1 priming', 'z-score priming decay']
        for i, ax in enumerate(axes):

            left = width + (bar_width*i) + (bar_pad*i)
            cbar_ax = fig.add_axes([left, bottom, bar_width, bar_height])
            cbar = plt.colorbar(ax, cax=cbar_ax)
            n = 5
            th = self.thresh[i]
            ticks = [round(th + th * (((s)/float(n-1))), 1) for s in range(n)]
            print ticks
            cbar.set_ticks(ticks)
            cbar.set_label(labels[i])
            plot_color(cbar_ax, color='#ffffff')

        fig.savefig(self.filename,
                    facecolor=fig.get_facecolor(),
                    edgecolor='none')


if __name__ == '__main__':

    root = '/mnt/data/FMRI/_testing/ForPublication/'
    model = 'decay_final2_flame12'
    con = 8
    fname = 'zstat1.nii.gz'
    contrasts = [os.path.join(root, model, str(con), '_fixedflameo0', fname),
                 os.path.join(root, model, str(con+1), '_fixedflameo0', fname)]
    cmaps = ['summer', 'winter']
    background = 'MNI152_T1_2mm.nii.gz'
    con2img = Con2Img('figure_6.png',
                      contrasts=contrasts,
                      fdr_corrected=True,
                      background=background,
                      cmaps=cmaps,
                      center=56,
                      masks=[os.path.join(root, 'mask.nii.gz'),
                             os.path.join(root, 'clr_mask.nii.gz')],
                      skip=1,
                      ntiles=12,
                      neg=True,
                      intersect=False,
                      thresh=2.3,
                      axis=2)

    con2img.save_image()
