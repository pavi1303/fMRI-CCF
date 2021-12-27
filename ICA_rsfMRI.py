import nilearn as nil
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
import os
import pandas as pd
from nilearn import datasets

#Importing datasets
rs_dataset = datasets.fetch_development_fmri(n_subjects=5)
func_filenames = rs_dataset.func

img_sample = nib.load('ds114_sub009_t2r1.nii')
n_voxels = np.prod(img_sample.shape[:-1])
n_trs = img_sample.shape[-1]
data = img_sample.get_data()
std_devs = []
for vol_no in range(n_trs):
   vol = data[..., vol_no]
   std_devs.append(np.std(vol))

plt.plot(std_devs)
voxels_by_time = data.reshape((n_voxels, n_trs))
std_devs_vectorized = np.std(voxels_by_time, axis=0)
assert np.allclose(std_devs_vectorized, std_devs)

#Trial 1 of importing all the docs
sub_id = str(8).zfill(3)
dir_path = 'Z:/COBRE_SCANS/008'
sub_path = 'sub_id'
path = os.path.join(dir_path,sub_path)
subdir = ''


for root, subdirs, files in os.walk('Z:/COBRE_SCANS/008'):
    for d in subdirs:
        if
pat_id = [8,13,23]
pat_id1 = str(pat_id)
pat_id_str = pat_id1.zfill(4)
subdir_list=[]
subdir_
for root, subdirectories, files in os.walk(dir_path):
    for subdirectory in subdirectories:
        print(os.path.join(root, subdirectory))

#Trial 2
basepath = 'Z:/COBRE_SCANS/008'
str_pat = 'cmrr_mbep2d_bold_AP_MB8_2mm_ISO'

for pattern in os.listdir(basepath):
    path = os.path.join(basepath,pattern)
    if os.path.isdir(path):
        if pattern == str_pat:
            print(path)