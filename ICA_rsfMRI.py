import os
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
import scipy as sp

# CENTERING BASED ON MEANS
def center_mean(x):
    mean = np.mean(x, axis=1, keepdims=True)
    centered = x - mean
    return centered
# COVARIANCE OF THE MATRIX
def cov(x):
    mean = np.mean(x, axis=1, keepdims=True)
    n = np.shape(x)[1] - 1
    m = x - mean
    return (m.dot(m.T)) / n
# WHITENING BASED ON SVD
def whiten(x):
    print('Performing whitening operation...')
    coVarM = cov(x)
    U, S, V = np.linalg.svd(coVarM)
    d = np.diag(1.0 / np.sqrt(S))
    whiteM = np.dot(U, np.dot(d, U.T))
    Xw = np.dot(whiteM, x)
    print('Whitening operation - DONE')
    return Xw, whiteM
# GET THE LIST OF NII FILES IN THE LOCATION
def _get_list_of_nii(location):
    if os.path.exists(location):
        os.chdir(location)
    else:
        print("Current working directory doesn't exist")
    list_of_nii=[]
    for files in os.listdir():
        if files.endswith(".nii"):
            list_of_nii.append(files)
    return list_of_nii
# GET THE DETAILS OF A TEMPLATE NII IMAGE
def _getnii_details(location,filename):
    if os.path.exists(location):
        os.chdir(location)
    else:
        print("Current working directory doesn't exist")
    img = nib.load(filename)
    data = img.get_fdata()
    n_voxels = np.prod(data.shape[:-1])
    n_trs = data.shape[-1]
    vols_shape = data.shape[:-1]
    nii_dim = data.shape
    aff = img.affine
    return nii_dim,vols_shape,n_voxels,n_trs,aff
# CONCATENATION OF TRS TO GET 4D NIFTI FILE
def _nii_concat(input_path, save_path):
    for dirpath,dirs, filenames in os.walk(input_path):
        for subdir in dirs:
            subpath = os.path.join(input_path, subdir)
            niilist = _get_list_of_nii(subpath)
            imgcat = nib.concat_images(niilist)
            nib.save(imgcat, os.path.join(save_path, 'MNI-' + subdir + '.nii'))
            del imgcat
# GENERATE 2D VOXEL-TIME SERIES DATA FROM 4D NIFTI IMAGE
def _obtain_vt_data(nii_image):
    n_voxels = np.prod(nii_image.shape[:-1])
    n_trs = nii_image.shape[-1]
    data = nii_image.get_fdata()
    vt_data = (data.reshape((n_voxels, n_trs)))
    return vt_data
# DO PCA & TEMPORAL CONCATENATION OF THE DATA
def _temporal_concat(location,n_comp,n_vxl):
    list_of_nii = _get_list_of_nii(location)
    length = len(list_of_nii)
    tempcat_dat = np.empty([1,n_vxl])
    for i in range(length):
        pat_img = nib.load(list_of_nii[i])
        voxtime_dat = _obtain_vt_data(pat_img).T
        print("Performing PCA using " + str(n_comp) + " components for subject " + str((i+1)) + "...")
        PCA_red = _do_PCA_v2(voxtime_dat,n_comp)
        print("Performing temporal concatenation of subject " + str((i+1)) + "...")
        tempcat_dat = np.append(tempcat_dat, PCA_red, axis=0)
        del pat_img,n_voxels,n_trs,data, voxtime_dat
    print('Temporal concatenation -- DONE')
    tempcat_dat = np.delete(tempcat_dat,0,0)
    return tempcat_dat
# PLOT TO DECIDE ON THE NUMBER OF PCA COMPONENTS TO USE
def _decide_PCA_comp (x):
    import numpy as np
    import matplotlib.pyplot as plt
    x_c = center_mean(x)
    Cov_x = cov(x_c)
    Ut, St, Vt = np.linalg.svd(Cov_x)
    cumsum = np.cumsum(St/np.sum(St))
    x_ax = range(1,len(cumsum)+1)
    plt.plot(x_ax,cumsum)
    plt.show()
    plt.xlabel("Number of components")
    plt.ylabel("Cumulative explained variance")
    plt.xticks(x_ax)
    plt.xlim([70,80])
# DO PCA ANALYSIS
def _do_PCA_v2(x,n_components):
    pca = PCA(n_components=n_components,whiten=False)
    pca_red = (pca.fit_transform(x.T)).T
    return pca_red
# OUTPUTTING THE RESULT OF ICA ANALYSIS AS NIFTI
def _save_ica_nifti(mat_loc,mat_filename,img_affine,vol_shape,dest_loc):
    os.chdir(mat_loc)
    ica_mat = loadmat(str(mat_filename) + '.mat')
    group_sm = ica_mat[mat_filename]
    # CONVERTING TO Z-SCORE MAPS
    group_sm -= group_sm.mean(axis=0)
    group_sm /= group_sm.std(axis=0)
    # THRESHOLDING IF NECESSARY
    #group_sm[np.abs(group_sm) > thresh] = 0
    # RESHAPING THE ICA MATRIX TO 4D
    group_sm_4D = group_sm.T.reshape(vol_shape+(group_sm.shape[0], ))
    for i in range(group_sm.shape[0]):
        gica_comp = group_sm_4D[...,i]
        gica_comp_img = nib.Nifti1Image(gica_comp,img_affine)
        if not os.path.exists(dest_loc):
            os.makedirs(dest_loc)
        nib.save(gica_comp_img, os.path.join(dest_loc, 'gICA_component_' + str(i + 1) + '.nii'))
        del gica_comp, gica_comp_img
# DUAL REGRESSION
def _dual_regression(location,gICA_sm):
    gICA_sm -= gICA_sm.mean(axis=0)
    gICA_sm /= gICA_sm.std(axis=0)
    sublist = _get_list_of_nii(location)
    number = len(sublist)
    for i in range(number):
        sub_img = nib.load(sublist[i])
        sub_vt = _obtain_vt_data(sub_img)
        sm = np.linalg.pinv(gICA_sm.T)
# END OF ALL THE FUNCTIONS FOR NOW


dim,vol,vox,trs,affi = _getnii_details(loc,file)
path = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
components = 20
pca_tcat = temporal_concat(path, components,vox)
savemat("pca_red.mat",{'pca_tcat':pca_tcat})


ica_red_mat = ica_mat[3]
_decide_PCA_comp(pca_tcat)
ICA = FastICA(whiten=False)
ica = ICA.fit_transform(pca_tcat.T).T
ICA_mat = _do_ICA(2,ICA_premat,10,z_score=True,thresh_method='max')

# IMPORTING THE MAT FILE THAT HAS THE SPATIAL MAPS
from scipy.io import loadmat
ica_mat = loadmat('ica_red.mat')
ICA_mat = ica_mat['X_reduced1']
# CONVERTING THESE SPATIAL MAPS TO Z-SCORE MAPS
# QN - When converting these into z-score maps, is it across different voxels within a component (or)
# across components per voxel
ICA_mat -= ICA_mat.mean(axis=0)
ICA_mat /= ICA_mat.std(axis=0)
# THRESHOLDING THESE ICA MAPS
ICA_mat[np.abs(ICA_mat) > 2] = 0
# RESHAPING THESE BACK TO 4D NIFTI FILES
i=0
ICA_MAT_4D = ICA_mat.T.reshape(vol+(ICA_mat.shape[0],))
for i in range(ICA_mat.shape[0]):
    ICA_img = ICA_MAT_4D[...,i]
    ICA_comp_img = nib.Nifti1Image(ICA_MAT_4D,affi)
    nib.


ni_img = nib.Nifti1Image(ICA_MAT_4D,affi)
nib.save(ni_img,os.path.join(path,'ICA_comp.nii'))
del ICA_img
# TRYING OUT DUAL REGRESSION
# STEP 1 - SPATIAL REGRESSION
# Using the template spatial maps as the
def _dual_regression(location,gICA_sm):
    gICA_sm -= gICA_sm.mean(axis=0)
    gICA_sm /= gICA_sm.std(axis=0)
    sublist = _get_list_of_nii(location)
    number = len(sublist)
    for i in range(number):
        sub_img = nib.load(sublist[i])
        sub_vt = _obtain_vt_data(sub_img)
        sm = np.linalg.pinv(gICA_sm.T)


ICA_mat -= ICA_mat.mean(axis=0)
ICA_mat /= ICA_mat.std(axis=0)


y = np.linalg.pinv(ICA_mat.T)
ts = np.dot(y,voxtime_dat)
ts = y*(voxtime_dat)

# Trying out my ica function

fileloc = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
name = 'ica_red_mat'
aff = img.affine
desloc = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/2.ICA/gICA_SM'

_save_ica_nifti(fileloc,name,aff,vol,desloc)

mat_loc = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
mat_filename = 'ica_red_mat'
img_affine = img.affine
dest_loc = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/2.ICA/gICA_SM'
vol_shape = vol
i=0




















































































































































