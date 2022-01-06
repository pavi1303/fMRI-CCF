import os
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from sklearn.decomposition import PCA
from sklearn.decomposition import FastICA
from scipy.io import loadmat,savemat

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

loc = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
file = 'MNI-008.nii'

dim,vol,vox,trs = _getnii_details(loc,file)

def _obtain_vt_data(nii_image):
    n_voxels = np.prod(nii_image.shape[:-1])
    n_trs = nii_image.shape[-1]
    data = nii_image.get_fdata()
    vt_data = (data.reshape((n_voxels, n_trs)))
    return vt_data

vt = _obtain_vt_data(img)



loc = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'

l = _get_list_of_nii(loc)


#--------------REFORMATTED FUNCTIONS--------------------#
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
# DO PCA ANALYSIS
def _do_PCA_v2(x,n_components):
    pca = PCA(n_components=n_components,whiten=True)
    pca_red = (pca.fit_transform(x.T)).T
    return pca_red
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


#----------------COPY OF FUNCTIONS FROM MY MAIN SCRIPT
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
    aff = img.affine
    return nii_dim,vols_shape,n_voxels,n_trs,aff
# CONCATENATION OF TRS TO GET 4D NIFTI FILE
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
# GENERATE THE TIME SERIES DATA
def _obtain_vt_data(nii_image):
    n_voxels = np.prod(nii_image.shape[:-1])
    n_trs = nii_image.shape[-1]
    data = nii_image.get_fdata()
    vt_data = (data.reshape((n_voxels, n_trs)))
    return vt_data
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
        voxtime_dat = _obtain_vt_data(pat_img).T
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

# DUAL REGRESSION
def _dual_regression(group_sm,img_affine,vol_shape,sub_loc,save_loc):
    group_sm -= group_sm.mean(axis=0)
    group_sm /= group_sm.std(axis=0)
    sublist = _get_list_of_nii(sub_loc)
    number = len(sublist)
    for i in range(number):
        ss_img = nib.load(sublist[i])
        ss_vt = _obtain_vt_data(ss_img)
        # Dual regression I - Spatial regression
        ss_tc = np.dot(np.linalg.pinv(group_sm.T), ss_vt)
        # Dual regression II - Temporal regression
        ss_sm = np.dot(sub_vt, np.linalg.pinv(ss_tc))
        # Conversion of these maps to z-score maps
        ss_sm -= ss_sm.mean(axis=0)
        ss_sm /= ss_sm.std(axis=0)
        subdir = str(sublist[i])
        ss_saveloc = os.path.join(save_loc,subdir)
        if not os.path.exists(ss_saveloc):
            os.makedirs(ss_saveloc)
        ss_sm_4D = sub_sm.T.reshape(vol_shape + (sub_sm.shape[0],))
        for j in range(sub_sm.shape[0]):
            ss_ica_comp = ss_sm_4D[..., j]
            ss_ica_comp_img = nib.Nifti1Image(ss_ica_comp, img_affine)
            nib.save(ss_ica_comp_img, os.path.join(ss_saveloc, 'ssICA_component_' + str(i + 1) + '.nii'))
            del ss_ica_comp, ss_ica_comp_img


# TRYING OUT SAVING EACH OF THE COMPONENTS OF ICA AS NIFTI
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
    nib.save(gica_comp_img,os.path.join(dest_loc, 'gICA_component_' + str(i+1) + '.nii'))
    del gica_comp,gica_comp_img
# TRYING MY CONCEPT OF CREATING DIRECTORIES
dest = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/2.ICA'
subdest = 'Group_spatial_maps'
final_dest = os.path.join(dest,subdest)
if not os.path.exists(final_dest):
    os.makedirs(final_dest)
