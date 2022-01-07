import os
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
import scipy as sp
from scipy.io import savemat
from scipy.io import loadmat

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
def _get_list_of_ext(location,ext):
    if os.path.exists(location):
        os.chdir(location)
    else:
        print("Current working directory doesn't exist")
    list_of_files = []
    for files in os.listdir():
        if files.endswith(ext):
            list_of_files.append(files)
    return list_of_files
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
            niilist = _get_list_of_ext(subpath, ".nii")
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
# DO SUBJECT WISE PCA AND SAVE THEM
def _subject_PCA(location,n_comp,save_location):
    list_of_nii = _get_list_of_ext(location,".nii")
    length = len(list_of_nii)
    for i in range(length):
        subloc = list_of_nii[i]
        pat_img = nib.load(subloc)
        voxtime_dat = _obtain_vt_data(pat_img).T
        print("Performing PCA using " + str(n_comp) + " components for subject " + str(subloc[4:7]) + "...")
        PCA_red = _do_PCA_v2(voxtime_dat,n_comp)
        if not os.path.exists(save_location):
            os.makedirs(save_location)
        np.save(os.path.join(save_location,'PCA_' + str(subloc[4:7] + '.npy')),PCA_red)
        print("PCA reduction done for subject " + str(subloc[4:7]) + ".")
        del pat_img, voxtime_dat, PCA_red
    print('PCA reduction -- DONE')
# PERFORM TEMPORAL CONCATENATION OF THE PCA COMPONENTS ACROSS SUBJECTS
def _temporal_concat(location,n_vxl):
    list_of_npy = _get_list_of_ext(location,".npy")
    length = len(list_of_npy)
    tempcat_dat = np.empty([1,n_vxl])
    for i in range(length):
        subloc = list_of_npy[i]
        pat_pca = np.load(subloc)
        print("Performing temporal concatenation for subject " + str(subloc[4:7]) + "...")
        tempcat_dat = np.append(tempcat_dat, pat_pca, axis=0)
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
    group_sm = sp.stats.zscore(group_sm,axis=1)
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
# DIMENSIONALITY REDUCTION OF THE ICA SPACE
def _shortlist_ICA(x,list_of_comp):
    x_red = x[list_of_comp,:]
    return x_red
# DUAL REGRESSION
def _dual_regression(group_sm,img_affine,vol_shape,sub_loc,save_loc):
    group_sm -= group_sm.mean(axis=0)
    group_sm /= group_sm.std(axis=0)
    group_sm1 = sp.stats.zscore(group_sm,axis=0)
    sublist = _get_list_of_ext(sub_loc,".nii")
    number = len(sublist)
    for i in range(number):
        ss_img = nib.load(sublist[i])
        ss_vt = _obtain_vt_data(ss_img)
        # Dual regression I - Spatial regression
        ss_tc = np.dot(np.linalg.pinv(group_sm.T), ss_vt)
        # Dual regression II - Temporal regression
        ss_sm = np.dot(ss_vt, np.linalg.pinv(ss_tc)).T
        # Conversion of these maps to z-score maps
        ss_sm = sp.stats.zscore(ss_sm,axis=1)
        subdir = str(sublist[i])
        subdir1 = subdir[4:7]
        ss_saveloc = os.path.join(save_loc,subdir1)
        if not os.path.exists(ss_saveloc):
            os.makedirs(ss_saveloc)
        ss_sm_4D = ss_sm.T.reshape(vol_shape + (ss_sm.shape[0],))
        for j in range(ss_sm.shape[0]):
            ss_ica_comp = ss_sm_4D[..., j]
            ss_ica_comp_img = nib.Nifti1Image(ss_ica_comp, img_affine)
            nib.save(ss_ica_comp_img, os.path.join(ss_saveloc, 'ssICA_component_' + str(j + 1) + '.nii'))
        del ss_ica_comp, ss_ica_comp_img
    del ss_img,ss_vt,ss_tc,ss_sm,subdir, subdir1, ss_saveloc

# END OF ALL THE FUNCTIONS FOR NOW
# FUNCTIONS TO CREATE
# 1. Minor tweak to the DR to save the individual TC's as well
# 2. Averaging the component maps across all the subjects
# Testing of my pipeline - dry run 1

# Getting the details of a template NIFTI image
path = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
file = 'MNI-008.nii'
dim, vol, vox, trs, aff = _getnii_details(path, file)
# Performing subject-wise PCA and temporal concatenation
components = 100
pca_tcat = temporal_concat(path, components,vox)
# Performing final stage of PCA on the concatenated matrix
pca_red_tcat = _do_PCA_v2(pca_tcat,100)
os.chdir('C:/Users/PATTIAP/Desktop/COBRE_VF/Results/1.PCA')
# Saving the result as a mat file for ICA analysis in MATLAB
savemat("pca_red.mat",{'pca_red_tcat':pca_red_tcat})
# Importing the ICA result from MATLAB and saving as NIFTI images
name = 'ica_red_mat'
desloc = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/2.ICA/gICA_SM'
_save_ica_nifti(path, name, aff, vol, desloc)
# Reducing the ICA dimensions space
# Need to create a function
# Performing dual regression
save_l = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/3.DR/Subject_spatialmaps'
_dual_regression(ICA_mat, aff, vol, path, save_l)




















































































































































