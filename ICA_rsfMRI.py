import os
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
from scipy.io import savemat
# CONCATENATION OF NII FILES
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
            nib.save(imgcat, os.path.join(save_path, 'MNI-' + subdir + '.nii'))
            del imgcat
#
# OBTAINING THE NII IMAGE DIMENSIONS FOR LATER USE
def _getnii_details(location,filename):
    import os
    import nibabel as nib
    import numpy as np
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
    return nii_dim,vols_shape,n_voxels,n_trs
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
# DO PCA ANALYSIS
def _do_PCA (x,n_components):
    # Centering the data
    x_cen = center_mean(x)
    # Calculation of covariance
    cov_x = cov(x_cen[0])
    # SVD decomposition to obtain eigenvalues and eigen vector
    U, S, V = np.linalg.svd(cov_x)
    # No need of sorting - SVD already orders S in descending order
    # Building the eigen vector matrix
    W = U[:,0:n_components]
    # Projection of the PCA components to the original data matrix
    pca = (np.dot(x.T,W)).T
    return pca
# DO TEMPORAL CONCATENATION OF THE DATA
def temporal_concat(location,n_comp,n_vxl):
    if os.path.exists(location):
        os.chdir(location)
    else:
        print("Current working directory doesn't exist")
    list_of_nii=[]
    for files in os.listdir():
        if files.endswith(".nii"):
            list_of_nii.append(files)
    length = len(list_of_nii)
    tempcat_dat = np.empty([1,n_vxl])
    for i in range(length):
        pat_img = nib.load(list_of_nii[i])
        n_voxels = np.prod(pat_img.shape[:-1])
        n_trs = pat_img.shape[-1]
        data = pat_img.get_fdata()
        voxtime_dat = (data.reshape((n_voxels, n_trs))).T
        print("Performing PCA using " + str(n_comp) + " components for subject " + str((i+1)) + "...")
        #PCA_red = _do_PCA(voxtime_dat,n_comp)
        PCA_bi = PCA(n_components=n_comp,whiten=True)
        PCA_red = (PCA_bi.fit_transform(voxtime_dat.T)).T
        print("Performing temporal concatenation of subject " + str((i+1)) + "...")
        tempcat_dat = np.append(tempcat_dat, PCA_red, axis=0)
        del pat_img,n_voxels,n_trs,data, voxtime_dat
    print('Temporal concatenation -- DONE')
    tempcat_dat = np.delete(tempcat_dat,0,0)
    return tempcat_dat
# DO GROUP ICA
def _do_ICA(threshold,x,n_comp,z_score = None,thresh_method='min'):
    from sklearn.decomposition import FastICA
    ICA = FastICA(n_components=n_comp,whiten=False, algorithm="parallel",max_iter=1000)
    ica = ICA.fit_transform(x.T).T
    if z_score == True:
        ica -= ica.mean(axis=0)
        ica /= ica.std(axis=0)
    else:
        pass
    if thresh_method == 'min':
        ica[np.abs(ica)<threshold]=0
    elif thresh_method == 'max':
        ica[np.abs(ica)>threshold]=0
    else:
        pass
    return ica
# Plotting function to decide on the number of principal components
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

dim,vol,vox,trs = _getnii_details(loc,file)
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
ICA_MAT_4D = ICA_mat.T.reshape(vol+(ICA_mat.shape[0],))
for i in range(ICA_mat.shape[0]):
    ICA_img = ICA_MAT_4D[:,:,:,i]
ICA_MAT_4D[:,:,:,i]
# TRYING OUT DUAL REGRESSION
# STEP 1 - SPATIAL REGRESSION
# Using the template spatial maps as the

































































































































































