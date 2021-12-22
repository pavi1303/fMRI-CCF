import nilearn as nil
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
import os
import pandas as pd
from nilearn import datasets
from nilearn import plotting

if os.path.exists("C:/Users/PATTIAP/Desktop/Dataset/MNI dataset"):
    os.chdir("C:/Users/PATTIAP/Desktop/Dataset/MNI dataset")
else:
    print("Current working directory doesn't exist")

#importing the dataset
newlist=[]
for files in os.listdir():
    if files.endswith(".nii"):
        newlist.append(files)
    else:
        continue
for i in range(5):
    newlist.pop(0)
img_cat = nil.image.concat_imgs(newlist)
nib.save(img_cat,'MNI-0008.nii')
#Reshaping 4D images to 2D arrays
img1 = nib.load('MNI-0008.nii')
img1.shape
num_voxels = np.prod(img1.shape[:-1])
n_timeseries = img1.shape[-1]
data = img1.get_data()
std_voxels = []
for tsr in range(n_timeseries):
    vol = data[..., tsr]
    std_voxels.append(np.std(vol))
voxels_by_time = data.reshape(std_voxels, data.shape[-1])
brain_mask = datasets.load_mni152_brain_mask()
plotting.plot_roi(brain_mask, cmap='Paired')



img = nib.load('MNI_012_00001.nii')
vol_shape = img.shape
img1 = nib.load('MNI_012_00002.nii')
n_voxels = np.prod(vol_shape)
current_dir_list =  os.listdir()
print(current_dir_list)