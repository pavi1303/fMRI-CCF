import os, shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import nibabel as nib
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import scipy as sp
from scipy.io import savemat
from scipy.io import loadmat
from nilearn.masking import apply_mask
from nilearn.masking import unmask

# ------------- List of sub functions ----------------- #
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
    list_of_files.sort()
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
# GENERATE 2D VOXEL-TIME SERIES DATA FROM 4D NIFTI IMAGE
def _obtain_vt_data(nii_image):
    n_voxels = np.prod(nii_image.shape[:-1])
    n_trs = nii_image.shape[-1]
    data = nii_image.get_fdata()
    vt_data = (data.reshape((n_voxels, n_trs)))
    return vt_data
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
# DIMENSIONALITY REDUCTION OF THE ICA SPACE
def _shortlist_ICA(x,list_of_comp):
    x_red = x[list_of_comp,:]
    return x_red
# -------------------- List of main functions --------------- #
# PERFORM NII CONCATENATION AND SINGLE SUBJECT PCA
def _mask_subject_PCA(input_path, save_path, mask_path, n_comp):
    os.chdir(mask_path)
    mask = nib.load('standard_binary.nii', mmap=False)
    for dirpath, dirs, filenames in os.walk(input_path):
        for subdir in dirs:
            subpath = os.path.join(input_path, subdir)
            os.chdir(subpath)
            niilist = _get_list_of_ext(subpath, ".nii")
            niilist.sort()
            niilist = niilist[15:]
            print('Performing nii concatenation for subject ' + str(subdir) + '...')
            pat_img = nib.concat_images(niilist)
            print('Applying mask and obtaining time series data for subject ' + str(subdir) + '...')
            os.chdir(mask_path)
            voxtime_dat = apply_mask(pat_img, mask)
            voxtime_dat = sp.stats.zscore(voxtime_dat, axis=0)
            print("Performing PCA using " + str(n_comp) + " components for subject " + str(subdir) + "...")
            PCA_red = _do_PCA_v2(voxtime_dat, n_comp)
            if not os.path.exists(save_path):
                os.makedirs(save_path)
            os.chdir(save_path)
            savemat('PCA_' + str(subdir) + '.mat', {'PCA_red': PCA_red})
            #np.save(os.path.join(save_path, 'PCA_' + str(subdir) + '.npy'), PCA_red)
            print("PCA reduction done for subject " + str(subdir) + ".")
            del pat_img, voxtime_dat, subpath, PCA_red
    print('PCA reduction -- DONE')
# PERFORM TEMPORAL CONCATENATION OF THE PCA COMPONENTS ACROSS SUBJECTS
def _temporal_concat(location,n_vxl, start, end):
    list_of_npy = _get_list_of_ext(location,".npy")
    list_of_npy.sort()
    list_of_npy = list_of_npy[start:end]
    length = len(list_of_npy)
    tempcat_dat = np.empty([1,n_vxl])
    for i in range(length):
        subloc = list_of_npy[i]
        pat_pca = np.load(subloc)
        #pat_pca = pat_pca.astype('float16')
        print("Performing temporal concatenation for subject " + str(subloc[4:7]) + "...")
        tempcat_dat = np.append(tempcat_dat, pat_pca, axis=0)
    print('Temporal concatenation -- DONE')
    tempcat_dat = np.delete(tempcat_dat,0,0)
    return tempcat_dat
# DO PCA ANALYSIS
def _do_PCA_v2(x,n_components):
    pca = PCA(n_components=n_components, whiten=False, svd_solver='full')
    pca_red = (pca.fit_transform(x.T)).T
    return pca_red
# OUTPUTTING THE RESULT OF ICA ANALYSIS AS NIFTI
def _save_ica_nifti(mat_loc,mat_filename,dest_loc, mask_path):
    os.chdir(mat_loc)
    ica_mat = loadmat(str(mat_filename) + '.mat')
    group_sm = ica_mat[mat_filename]
    # CONVERTING TO Z-SCORE MAPS
    group_sm = sp.stats.zscore(group_sm, axis=0)
    # RESHAPING THE ICA MATRIX TO 4D
    os.chdir(mask_path)
    mask = nib.load('standard_binary.nii', mmap=False)
    group_sm_4D = unmask(group_sm, mask)
    group_sm_data = group_sm_4D.get_fdata()
    #group_sm_4D = group_sm.T.reshape(vol_shape+(group_sm.shape[0], ))
    for i in range(group_sm.shape[0]):
        gica_comp = group_sm_data[...,i]
        gica_comp_img = nib.Nifti1Image(gica_comp, mask.affine)
        if not os.path.exists(dest_loc):
            os.makedirs(dest_loc)
        nib.save(gica_comp_img, os.path.join(dest_loc, 'gICA_component_' + str(i + 1) + '.nii'))
        del gica_comp, gica_comp_img
# DUAL REGRESSION
def _dual_regression(tc_loc, sm_loc, sub_path, mat_loc, mat_filename):
    os.chdir(mat_loc)
    ica_mat = loadmat(str(mat_filename) + '.mat')
    group_sm = ica_mat[mat_filename]
    group_sm = sp.stats.zscore(group_sm, axis=1)
    os.chdir(mask_path)
    mask = nib.load('MNI_152_mask_v2.nii', mmap=False)
    for dirpath, dirs, filenames in os.walk(sub_path):
        for subdir in dirs:
            subpath = os.path.join(sub_path, subdir)
            niilist = _get_list_of_ext(subpath, ".nii")
            niilist.sort()
            niilist = niilist[15:]
            print('Performing nii concatenation for subject ' + str(subdir) + '...')
            pat_img = nib.concat_images(niilist)
            print('Applying mask and obtaining time series data for subject ' + str(subdir) + '...')
            ss_vt = apply_mask(pat_img, mask)
            # Dual regression I - Spatial regression
            print('Doing spatial regression for subject ' + str(subdir) + '...')
            ss_tc = np.dot(np.linalg.pinv(group_sm.T), ss_vt)
            if not os.path.exists(tc_loc):
                os.makedirs(tc_loc)
            np.save(os.path.join(tc_loc, 'ss_tc_' + str(subdir) + '.npy'), ss_tc)
            print('Spatial regression -- DONE for subject ' + str(subdir) + '.')
            # Dual regression II - Temporal regression
            print('Doing temporal regression for subject ' + str(subdir) + '...')
            ss_sm = np.dot(ss_vt, np.linalg.pinv(ss_tc)).T
            # Conversion of these maps to z-score maps
            ss_sm = sp.stats.zscore(ss_sm, axis=1)
            ss_saveloc = os.path.join(sm_loc, subdir)
            if not os.path.exists(ss_saveloc):
                os.makedirs(ss_saveloc)
            ss_sm_4D = unmask(ss_sm, mask)
            for j in range(ss_sm.shape[0]):
                ss_ica_comp = ss_sm_4D[..., j]
                ss_ica_comp_img = nib.Nifti1Image(ss_ica_comp, mask.affine)
                nib.save(ss_ica_comp_img, os.path.join(ss_saveloc, 'ssICA_component_' + str(j + 1) + '.nii'))
                del ss_ica_comp, ss_ica_comp_img
            print('Temporal regression -- DONE for subject ' + str(subdir) + '.')
            del ss_vt, ss_tc, ss_sm, pat_img, ss_saveloc
    print('Dual regression -- DONE')


# END OF ALL THE FUNCTIONS FOR NOW
# FUNCTIONS TO CREATE
# 1. Function to generate the correlation matrices. See if we can include it in the dual regression code.
# Generate a separate function for calculating FC matrix subject wise to yield an M X M matrix (M-ICA components)
# 2. Averaging the component maps across all the subjects
# 3. How to calculate partial correlation coefficient
#----------------------------------------------------------------------#
# Modifications to be made to the dual regression function
# (b) Calculate the Pearson Correlation Coefficient save the array, images
# (c) See if you can calculate the Partial Correlation Coefficient in the same way
#----------------------------------------------------------------------#
# Parameters for the ICA function in MATLAB
#method=11; a1=1; var_normal=1; eps=1E-6; A0=[]; shift=[]; Sigma=0; determine_flip=1;

#----------------------My run 1-----------------------#
# Path variables for all the save locations
sub_loc = 'E:/LRCBH/Data/COBRE-MNI/Trial'
mask_loc = 'E:/LRCBH/MNI_segmented'
ss_pca = 'E:/LRCBH/Results/1.PCA/With_standard/Trial'
ica_result = 'E:/LRCBH/Results/2.ICA/With_standard/v2'
temp_nii = 'E:/LRCBH/COBRE-MNI'
pca_result = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/1.PCA/Final'
ss_tc = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/3.DR/Subject_timecourses'
ss_sm = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/3.DR/Subject_spatialmaps'
#--------------------COPYING NII FILES-----------------------------#

#------------------------PCA REDUCTION------------------------------#
# Performing subject wise dimensionality reduction - PCA
comp = 200
_mask_subject_PCA(sub_loc, ss_pca, mask_loc, comp)

# Importing the ICA result from MATLAB and saving as NIFTI images
name = 'gICA_30_v2'
_save_ica_nifti(ica_result, name, ica_result, mask_loc)

# Getting the details from a template NIFTI image
file = 'MNI-008.nii'
dim, vol, vox, trs, aff = _getnii_details(temp_nii, file)
# Performing temporal concatenation on the obtained subject PCA's
vox = 187997
pca_tcat5 = _temporal_concat(ss_pca, vox, 8, 10)
os.chdir(ss_pca)
savemat('pca_tcat_trial5.mat', {'pca_tcat5': pca_tcat5})
#np.save(os.path.join(pca_result, 'pca_tcat2.npy'), pca_tcat)
#Saving the matlab variables in parts

# Need to save the tcat variables into different parts of .mat variables
# Performing final stage of PCA on the concatenated matrix - if necessary
#Saving both the python and the matlab versions of the result
#pca_red_tcat = _do_PCA_v2(pca_tcat, 100)
#np.save(os.path.join(final_pca, 'pca_red_tcat.npy'), pca_red_tcat)
#os.chdir(final_pca)
#savemat('pca_red_tcat.mat', {'pca_red_tcat': pca_red_tcat})

# Performing dual regression
_dual_regression(aff, vol, ss_tc, ss_sm, sub_loc, ica_result, name)

###############My final ray of hope###################

# Redundant code - copying the files
root_folder = 'Z:/COBRE_SCANS'
smooth_folder = 'intermediate'
dest_root = 'E:/LRCBH/COBRE-MNI/Individual_data'
os.chdir('C:/Users/PATTIAP/Documents')
df = pd.read_excel('Bold_location.xlsx', sheet_name='Bold')
bold_folders = df['Names'].tolist()
pat_loc = []
for dirpath, dirs, filenames in os.walk(dest_root):
    pat_loc.append(dirs)
    pat_folders = pat_loc[0]
    pat_folders.sort()
length = len(pat_folders)
# Performing the copying operation
for i in range(length):
    input_path = os.path.join(root_folder, pat_folders[i], bold_folders[i], smooth_folder)
    save_path = os.path.join(dest_root, pat_folders[i])
    niilist = []
    if not os.path.exists(input_path):
        print("Input path doesnt exist")
    if not os.path.exists(save_path):
        print("Input path doesnt exist")
    for files in os.listdir(input_path):
        if files.endswith('.nii') and files.startswith('s'):
            niilist.append(files)
            niilist.sort()
    os.chdir(input_path)
    print('Copying the nii files for subject ' + str(pat_folders[i] + '...'))
    for f in niilist:
        shutil.copy(f,save_path)
    del niilist, input_path, save_path
print('Copying operation -- DONE')


#CHECKING MY PCA OPERATION
os.chdir('E:/LRCBH/Old/Results/1.PCA/Trial')
x = loadmat('vt_008.mat')
x_mean = np.mean(data, axis =0)
data = x['voxtime_dat']
c = np.cov(data)
u, s, v = np.linalg.svd(c)
W = u[:, 0:200]
# Projection of the PCA components to the original data matrix
pca = (np.dot(data.T, W)).T
data_red = _do_PCA_v2(data, 200)

os.chdir('E:/ICA_test_sub008_2017')
Ydata = loadmat('Ydata.mat')
dat = Ydata['Ydata']
c = np.cov(dat)
u, s, v = np.linalg.svd(c)
W = u[:, 0:100]
# Projection of the PCA components to the original data matrix
pca = (np.dot(dat.T, W)).T

# Creating a group scatter plot with regression line 



































































































































