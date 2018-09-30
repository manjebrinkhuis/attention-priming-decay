from nipype.interfaces import fsl, freesurfer
import numpy as np
import subprocess as sp
import nibabel as nib
import os
import matplotlib.pyplot as plt
from scipy import stats


def dict2opts(dct, doubledash=True,
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
            opts.append('{0}{1}{2}{3}'.format(dash, opt, assign_sign, dct[kw]))
        elif dct[kw] is True:
            opts.append('{0}{1}'.format(dash, opt))

    return opts


def compare_coords(coord1, coord2, dist):
    """
    Return boolean. True if distance between
    2 coords is smaller than dist.
    """
    x, y, z = np.array(coord1) - np.array(coord2)
    return (x**2+y**2+z**2)**.5 <= dist


def group_coords(coords, min_dist, **kwargs):

    """
    Return list of lists. Recursively compare coords and
    return groups of coords that are within range of min_dist
    of eachother.
    """

    if 'groups' not in kwargs.keys():
        groups = []
    else:
        groups = kwargs['groups']

    if 'indices' in kwargs:
        indices = kwargs['indices']
    else:
        indices = range(len(coords))

    if coords:

        compare = coords[0]
        rest = coords[1:]

        index = indices[0]
        rest_indices = indices[1:]
        groups.append([(index, compare)])

        include = []
        for i, coord in enumerate(rest):
            if compare_coords(compare, coord, min_dist):
                index = rest_indices[i]
                groups[-1].append((index, coord))

            else:
                include.append(i)

        rest = [rest[i] for i in include]
        indices = [rest_indices[i] for i in include]

        return group_coords(rest, min_dist,
                            groups=groups,
                            indices=indices)
    else:
        return groups


def get_radii(coords):
    """ Radii of x,y,z arrays in array. Distance from (0,0,0). """
    return [(x**2+y**2+z**2)**.5 for x, y, z in coords]


def make_sphere(radius=3):
    """
    Create 3D array with ones and zeros. 1 if voxel
    falls within sphere with given radius, else 0.
    """
    box = np.ones((radius*2+1,)*3)
    d = radius*2+1
    xyzs = [(x, y, z) for x in range(d) for y in range(d) for z in range(d)]
    coords = np.array(xyzs) - radius
    dist = get_radii(coords)

    for i, (x, y, z) in enumerate(xyzs):

        if dist[i] > (radius):
            box[x, y, z] = 0
        else:
            box[x, y, z] = 1

    return box


class FSLCluster():
    """ Use FSL to create clusters. """

    def __init__(self, in_file, out_file=None, thresh=1.96, **kwargs):
        self.in_file = in_file
        self.out_file = out_file
        self.thresh = thresh
        self.cmd = ['cluster', '--in='+in_file,
                    '--thresh='+str(thresh)]

        self.cmd += dict2opts(kwargs)

    def run(self):

        # Run subprocess
        pipe = sp.Popen(' '.join(self.cmd),
                        shell=True,
                        stdout=sp.PIPE)

        # Process output
        lines = pipe.stdout.read().split('\n')

        # Analyze header
        hdr = lines[0].split('\t')

        # Create tuples for data recarray
        data = [tuple(line.split('\t')) for line in lines[1:-1]]

        if self.out_file:
            np.savetxt(self.out_file, data, delimiter='\t', fmt='%s', header='\t'.join(hdr))

        return np.array(data, dtype=[(name, 'float') for name in hdr])


class Cluster():
    def __init__(self, source, target, min_dist=2.5, sphere_radius=4, **kwargs):
        self.source = source
        self.target = target

        path, fname = os.path.split(self.source)
        self.index_file = os.path.join(path, 'out_index.nii.gz')

        fsl_cluster = FSLCluster(source, oindex=self.index_file, **kwargs)
        self.clusters = fsl_cluster.run()

        self.min_dist = min_dist
        self.sphere_radius = sphere_radius
        self.set_groups(min_dist)

    def set_groups(self, min_dist):
        """ """
        self.min_dist = min_dist
        coords = self.get_coords()

        groups = group_coords(coords, min_dist)
        indices = [[coord[0] for coord in group] for group in groups]

        self.groups = indices

    def get_coords(self, cluster_index=None, pkt='MAX'):
        """ """
        fields = [pkt + ' X (vox)', pkt + ' Y (vox)', pkt + ' Z (vox)']
        coords = [list(coord) for coord in self.clusters[fields]]
        indices = cluster_index if not cluster_index is None else range(len(coords))
        return [coords[i] for i in indices]

    def get_spheres(self):

        zstat1 = nib.load(self.source)
        rad = self.sphere_radius
        sphere = make_sphere(rad)

        path, fname = os.path.split(self.source)
        masks = []

        for i, group in enumerate(self.groups):

            mask = np.zeros(zstat1.shape, dtype=bool)

            for index in group:

                coord = self.get_coords()[index]
                x, y, z = coord
                x_slice = slice(x-rad, x+rad+1)
                y_slice = slice(y-rad, y+rad+1)
                z_slice = slice(z-rad, z+rad+1)
                mask[x_slice, y_slice, z_slice] += sphere.astype(bool)

            mask_fname = os.path.join(path, 'mask_'+str(i)+'.nii.gz')
            mask_img = nib.Nifti1Image(mask.astype(int),
                                       affine=zstat1.affine)
            nib.save(mask_img, mask_fname)
            masks.append(mask_fname)

        return masks

    def get_average_at_spheres(self, cluster_index=None, frac=.1):

        tgt_data = nib.load(self.target).get_data()
        avgs = []
        masks = self.get_spheres()
        for m in masks:
            mask = nib.load(m).get_data().astype(bool)
            num = len(tgt_data[mask])
            data = sorted(tgt_data[mask][:int(num*frac)])

            print int(num*frac), '/', num
            print np.mean(data), np.mean(tgt_data[mask])

            avg = np.mean(data)
            avgs.append(avg)

        return avgs

    def get_volume(self, cluster_index=None):
        """ """
        indices = cluster_index if cluster_index else range(len(self.clusters))
        return self.clusters['Voxels'][indices]

    def get_cluster_info(self, cluster_index=None):
        """ """
        indices = cluster_index if cluster_index else range(len(self.clusters))
        return [self.clusters[i] for i in indices]

    def get_value_at_cluster_peak(self, cluster_index=None):
        coords = self.get_coords(cluster_index=None, pkt='MAX')
        tgt_data = nib.load(self.target).get_data()
        values = []
        for coord in coords:
            x, y, z = coord
            values.append(tgt_data[x, y, z])

        return values

    def get_peak_voxel(self, cluster_index=None, _min=True):
        """ """

        tgt_data = nib.load(self.target).get_data()
        mask_data = nib.load(self.index_file).get_data()

        if cluster_index:
            indices = [cluster_index]
        else:
            indices = range(1, np.max(mask_data)+1)

        peaks = []
        for index in reversed(indices):
            mask = mask_data == index
            data = tgt_data[mask]
            peak = np.min(data) if _min else np.max(data)
            peaks.append(peak)

        return peaks

    def get_average_voxel(self, cluster_index=None):

        tgt_data = nib.load(self.target).get_data()
        mask_data = nib.load(self.index_file).get_data()

        if cluster_index:
            indices = [cluster_index]
        else:
            indices = range(1, np.max(mask_data)+1)

        avgs = []
        for index in reversed(indices):
            mask = mask_data == index
            data = tgt_data[mask]
            avg = np.mean(data)
            avgs.append(avg)

        return avgs


contrasts = [4,     # loc_lag12
             5,     # loc_lag1
             6,     # clr_lag2
             7,     # clr_lag12
             8,     # clr_lag1
             9,     # clr_lag2
             10]     # stim onset

path = os.path.join('/mnt/data/FMRI/_testing/ForPublication',
                    'decay_final2_flame12',
                    '{con}', '_fixedflameo0', 'zstat1.nii.gz')
out_fname = '{fx}_decay_con_{con}.nii.gz'
stats_fname = '{fx}_decay_con_{con}.txt'
thresh = 2.33
min_size = 40

reg_file = '/mnt/data/Apps/Linux/freesurfer/average/mni152.register.dat'

for con in contrasts:

    fname = path.format(con=con)
    cluster_min = FSLCluster(in_file=fname,
                             thresh=-thresh,
                             out_file=stats_fname.format(con=con, fx='min'),
                             _min=True,
                             minextent=min_size,
                             oindex=out_fname.format(con=con, fx='min'))
    cluster_min.run()

    cluster_max = FSLCluster(in_file=fname,
                             thresh=thresh,
                             out_file=stats_fname.format(con=con, fx='max'),
                             _min=False,
                             minextent=min_size,
                             oindex=out_fname.format(con=con, fx='max'))
    cluster_max.run()

    out_surf = '/mnt/data/FMRI/_testing/ForPublication/'\
        '{hemi}.decay_con_{con}.mgz'

    hemi = 'lh'
    fwhm = 1
    vol2surf = ['mri_vol2surf',
                '--mov', fname,
                '--surf-fwhm', str(fwhm),
                '--hemi', hemi,
                '--mni152reg',
                '--o', out_surf.format(con=con, hemi=hemi)]
    # sp.call(vol2surf)

    hemi = 'rh'
    vol2surf = ['mri_vol2surf',
                '--mov', fname,
                '--surf-fwhm', str(fwhm),
                '--hemi', hemi,
                '--mni152reg',
                '--o', out_surf.format(con=con, hemi=hemi)]
    # sp.call(vol2surf)
