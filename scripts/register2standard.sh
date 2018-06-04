struct="/data/ext/PhD/Studies/PrimingDecay/data/in_analysis/fs_subjects/sub_001/mri/orig.nii.gz"
brain="/data/ext/PhD/Studies/PrimingDecay/data/in_analysis/fs_subjects/sub_001/mri/brain.nii.gz"
MNI="/data/ext/PhD/Studies/PrimingDecay/data/in_analysis/nii/standard/MNI152_T1_2mm_brain.nii.gz"
omat2="aff_struct2mni.mat"
warp_struct2mni="warp_struct2mni.nii.gz"

SUBJECTS_DIR="/data/ext/PhD/Studies/PrimingDecay/data/in_analysis/fs_subjects"
inplane="/data/ext/PhD/Studies/PrimingDecay/data/in_analysis/nii/sub_001/ses_000/anatomy/inplane.nii.gz"
omat1="inplane2struct.mat"
bbreg="inplane_bbreg.dat"
betted="betted.nii.gz"
bet ${inplane} ${betted}
bbregister \
  --t1 \
  --6 \
  --fsl-dof 6 \
  --init-fsl \
  --fslmat ${omat1} \
  --reg ${bbreg} \
  --mov ${betted} \
  --o bbreg_out.nii.gz \
  --s sub_001

# tkregisterfv \
#     --mov /data/ext/PhD/Studies/PrimingDecay/data/in_analysis/nii/sub_001/ses_000/anatomy/inplane.nii.gz \
#     --reg inplane_bbreg.dat \
#     --surfs

tkregister2 --mov ${inplane} \
  --reg ${bbreg}
  --surf

# flirt \
#     -ref ${MNI} \
#     -in ${brain}  \
#     -omat ${omat2}
#
# fnirt \
#     --ref=${MNI} \
#     --in=${struct} \
#     --aff=${omat2} \
#     --cout=${warp_struct2mni}
#
# applywarp \
#     --ref=${MNI} \
#     --in=${inplane} \
#     --out=inplaneMNI.nii.gz \
#     --warp=${warp_struct2mni} \
#     --premat=${omat1}
