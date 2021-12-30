import os
import numpy as np
import nibabel as nib

#---------------------------FUNCTIONS--------------------------------------#
#------------------Import and concat the nii files-------------------------#
def nii_concat(input_path, save_path):
    for dirpath,dirs, filenames in os.walk(input_path):
        for subdir in dirs:
            niilist = []
            subpath = os.path.join(input_path, subdir)
            if os.path.exists(subpath):
                os.chdir(subpath)
            else:
                print("Current working directory doesn't exist")
            for files in os.listdir():
                if files.endswith(".nii"):
                    niilist.append(files)
            imgcat = nib.concat_images(niilist)
            nib.save(imgcat,os.path.join(save_path,'MNI-'+subdir+'.nii'))
            del imgcat
#-------------------Temporal concatenation of the data--------------------#
def temporal_concat(location):
    if os.path.exists(location):
        os.chdir(location)
    else:
        print("Current working directory doesn't exist")
    list_of_nii=[]
    for files in os.listdir():
        if files.endswith(".nii"):
            list_of_nii.append(files)
    length = len(list_of_nii)
    tempcat_dat = np.empty([1,902629])
    for i in range(length):
        pat_img = nib.load(list_of_nii[i])
        n_voxels = np.prod(pat_img.shape[:-1])
        n_trs = pat_img.shape[-1]
        data = pat_img.get_fdata()
        voxtime_dat = (data.reshape((n_voxels, n_trs))).transpose()
        tempcat_dat = np.append(tempcat_dat, voxtime_dat, axis=0)
        del pat_img,n_voxels,n_trs,data, voxtime_dat
    tempcat_dat = np.delete(tempcat_dat,0,0)
    return tempcat_dat
#------------------Preprocessing for multi-subject ICA-------------------#
def center_mean(x):
    mean = np.mean(x, axis=0, keepdims=True)
    centered = x - mean
    return centered, mean
# 2. Whitening operation based on covariance matrix
def cov(x):
    mean = np.mean(x, axis=0, keepdims=True)
    n = np.shape(x)[1] - 1
    m = x - mean
    return (m.dot(m.T)) / n
def whiten(x):
    coVarM = cov(x)
    U, S, V = np.linalg.svd(coVarM)
    d = np.diag(1.0 / np.sqrt(S))
    whiteM = np.dot(U, np.dot(d, U.T))
    Xw = np.dot(whiteM, x)
    return Xw, whiteM

source = 'C:/Users/PATTIAP/Desktop/Dataset/MNI dataset'
destination = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI/mni'

nii_concat(source,destination)
ICA_mat = temporal_concat(destination)

os.getcwd()
