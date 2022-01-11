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
        ss_tc = np.dot(np.linalg.pinv(group_sm.T), ss_vt.T)
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
def _dual_regression(group_sm,img_affine,vol_shape,sub_loc,save_loc):
    group_sm -= group_sm.mean(axis=0)
    group_sm /= group_sm.std(axis=0)
    group_sm1 = sp.stats.zscore(group_sm,axis=0)
    sublist = _get_list_of_nii(sub_loc)
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
        np.seterr(invalid='ignore')
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
# TRYING MY CONCEPT OF CREATING DIRECTORIES
dest = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/2.ICA'
subdest = 'Group_spatial_maps'
final_dest = os.path.join(dest,subdest)
if not os.path.exists(final_dest):
    os.makedirs(final_dest)

s = 'MNI-113.nii'

s1 = str(s[4:7])

group_sm = ICA_mat
img_affine = affi
sub_loc = l
dest_loc = save_l
i=0

ica_red_mat = ica_mat[3]
_decide_PCA_comp(pca_tcat)
ICA = FastICA(whiten=False)
ica = ICA.fit_transform(pca_tcat.T).T
ICA_mat = _do_ICA(2,ICA_premat,10,z_score=True,thresh_method='max')

# IMPORTING THE MAT FILE THAT HAS THE SPATIAL MAPS
from scipy.io import loadmat
os.chdir('C:/Users/PATTIAP/Desktop/COBRE_VF/Results/2.ICA')
ica_mat = loadmat('ica_red_mat.mat')
ICA_mat = ica_mat['ica_red_mat']
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
# Trying out my dual regression function

l = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
vol_shape = vol
save_l = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/3.DR/Subject_spatialmaps'

_dual_regression(ICA_mat,affi,vol,l,save_l)

mat_loc = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
mat_filename = 'ica_red_mat'
img_affine = img.affine
dest_loc = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/2.ICA/gICA_SM'
vol_shape = vol
i=0

list = [1,2,3,4,5,6]

ICA_mat_red = ICA_mat[list,:]

# Breakup the temporal concat function into two parts
# Generate a mini function for ICA space reduction
# Minor tweak to the dual regression function to save the tc's as well
def _shortlist_ICA(x,list_of_comp):
    x_red = x[list_of_comp,:]
    return x_red
ICA_x = _shortlist_ICA(ICA_mat,list)
def _get_list_of_ext(location, ext):
        if os.path.exists(location):
            os.chdir(location)
        else:
            print("Current working directory doesn't exist")
        list_of_files = []
        for files in os.listdir():
            if files.endswith(ext):
                list_of_files.append(files)
        return list_of_files
q = _get_list_of_ext(path,".mat")
# Breaking the temporal concat into two functions
def _temporal_concat(location,n_comp,n_vxl):
    list_of_nii = _get_list_of_ext(location,".nii")
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

comp = 10
s_l = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/1.PCA'
_subject_PCA(path, comp, s_l)
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

tcat = _temporal_concat(s_l,vox)
input_loc = 'E:/LRCBH/COBRE-MNI/Trial'
save_loc = 'E:/LRCBH/Results/1.PCA'
comp = 5

_concat_subject_PCA(comp, input_loc, save_loc)

# My new function for dual regression
os.chdir('D:/LRCBH/Results/2.ICA')
mat = loadmat('ica_red_mat.mat')
mat_file = mat['ica_red_mat']
tc_path='D:/LRCBH/Results/3.DR/Variables/Subject_timecourses'
sm_path='D:/LRCBH/Results/3.DR/Variables/Subject_spatialmaps'
sub_path='D:/LRCBH/COBRE-MNI/Trial'
def _dual_regression(img_affine, vol_shape, tc_loc, sm_loc, sub_path, mat_loc, mat_filename):
    os.chdir(mat_loc)
    ica_mat = loadmat(str(mat_filename) + '.mat')
    group_sm = ica_mat[mat_filename]
    group_sm = sp.stats.zscore(group_sm,axis=1)
    for dirpath, dirs, filenames in os.walk(sub_path):
        for subdir in dirs:
            subpath = os.path.join(sub_path, subdir)
            niilist = _get_list_of_ext(subpath, ".nii")
            print('Performing nii concatenation for subject ' + str(subdir) + '...')
            pat_img = nib.concat_images(niilist)
            print('Generating voxel time series data for subject ' + str(subdir) + '.')
            ss_vt = _obtain_vt_data(pat_img)
            # Dual regression I - Spatial regression
            print('Doing spatial regression for subject ' + str(subdir) + '...')
            ss_tc = np.dot (np.linalg.pinv (group_sm.T), ss_vt)
            if not os.path.exists(tc_loc):
                os.makedirs(tc_loc)
            np.save(os.path.join(tc_loc,'ss_tc_' + str(subdir) + '.npy'), ss_tc)
            print('Spatial regression -- DONE for subject ' + str(subdir) + '.')
            # Dual regression II - Temporal regression
            print('Doing temporal regression for subject ' + str (subdir) + '...')
            ss_sm = np.dot (ss_vt, np.linalg.pinv (ss_tc)).T
            # Conversion of these maps to z-score maps
            ss_sm = sp.stats.zscore (ss_sm, axis=1)
            ss_saveloc = os.path.join (sm_loc, subdir)
            if not os.path.exists(ss_saveloc):
                os.makedirs(ss_saveloc)
            ss_sm_4D = ss_sm.T.reshape (vol_shape + (ss_sm.shape[0],))
            for j in range (ss_sm.shape[0]):
                ss_ica_comp = ss_sm_4D[..., j]
                ss_ica_comp_img = nib.Nifti1Image (ss_ica_comp, img_affine)
                nib.save (ss_ica_comp_img, os.path.join (ss_saveloc, 'ssICA_component_' + str (j + 1) + '.nii'))
                del ss_ica_comp, ss_ica_comp_img
            print('Temporal regression -- DONE for subject ' + str(subdir) + '.')
            del ss_vt,ss_tc,ss_sm,pat_img,ss_saveloc


_dual_regression(mat_file,aff,vol,tc_path,sm_path,sub_path)

def _fullcorr(input_path, save_path):
    list_of_npy = _get_list_of_ext(input_path, ".npy")
    list_of_npy.sort()
    length = len(list_of_npy)
    for i in range(length):
        sub_tc = list_of_npy[i]
        pat_tc = np.load(sub_tc)
        corr = np.corrcoef(pat_tc)
        cm = sns.heatmap(corr, annot=True, vmin=-1, vmax=1, annot_kws={"size": 7}, cmap="YlGnBu")
        fig = cm.get_figure()
        os.chdir(save_path)
        fig.savefig(('FC_matrix_' + str(sub_tc[6:9]) + '.png'),dpi=300)
        del sub_tc, pat_tc, corr, cm, fig

x = str('ss_tc_008')
_fullcorr(ss_tc, ss_tc)

save_path = 'D:/LRCBH/Results/1.PCA/v2'
np.save(os.path.join(save_path, 'pca_red_tcat.npy'), pca_red_tcat)
arr = np.random.rand(10,850)
c = np.corrcoef(arr)
import seaborn as sns
plt.use("TkAgg")
del hm,fig
hm2 = sns.heatmap(c,annot=True,vmin=-1,vmax=1,annot_kws={"size": 7},cmap="YlGnBu")
fig2 = hm.get_figure()
fig2.savefig('Corrmat1.png',dpi=300)
plt.show()

def _gen_fullcorr(sub_tc, corr_png, corr_npy):
    mat_map = sns.heatmap(c,annot=True,vmin=-1,vmax=1,annot_kws={"size": 7})

os.chdir('E:/LRCBH/Results/1.PCA/v2')
x = np.load('pca_tcat.npy')
x1 = x[900:,:]
pca_tcat4 = x1
pca_tcat = pca_tcat.astype('float16')
savemat('pca_tcat_part4.mat',{'pca_tcat4': pca_tcat4})


# Dump from main script - Dated 1/10/2022 11:27 AM
# Performing nii concatenation of all the subjects
input_loc = 'E:/LRCBH/COBRE-MNI/Grp1-Controls'
save_loc = 'E:/LRCBH/Concatenated/Grp1-controls'
_nii_concat(input_loc, save_loc)
# Getting the details of a template NIFTI image
path = 'E:/LRCBH/Concatenated'
file = 'MNI-008.nii'
dim, vol, vox, trs, aff = _getnii_details(path, file)
# Performing subject-wise PCA and temporal concatenation
components = 100
pca_tcat = temporal_concat(path, components, vox)
# Performing final stage of PCA on the concatenated matrix
pca_red_tcat = _do_PCA_v2(pca_tcat, 100)
os.chdir('C:/Users/PATTIAP/Desktop/COBRE_VF/Results/1.PCA')
# Saving the result as a mat file for ICA analysis in MATLAB
savemat("pca_red.mat", {'pca_red_tcat': pca_red_tcat})
# Importing the ICA result from MATLAB and saving as NIFTI images
name = 'ica_red_mat'
desloc = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/2.ICA/gICA_SM'
_save_ica_nifti(path, name, aff, vol, desloc)
# Reducing the ICA dimensions space
# Need to create a function
# Performing dual regression
save_l = 'C:/Users/PATTIAP/Desktop/COBRE_VF/Results/3.DR/Subject_spatialmaps'
_dual_regression(ICA_mat, aff, vol, path, save_l)
# FUNCTIONS NOT IN USE
# CONCATENATION OF TRS TO GET 4D NIFTI FILE
def _nii_concat(input_path, save_path):
    for dirpath, dirs, filenames in os.walk(input_path):
        for subdir in dirs:
            subpath = os.path.join(input_path, subdir)
            niilist = _get_list_of_ext(subpath, ".nii")
            print('Performing nii concatenation for subject ' + str(subdir) + '...')
            imgcat = nib.concat_images(niilist)
            nib.save(imgcat, os.path.join(save_path, 'MNI-' + str(subdir) + '.nii'))
            del imgcat
            print('nii concatenation done for subject ' + str(subdir) + '.')
    print('NII concatenation -- DONE')
# DO SUBJECT WISE PCA AND SAVE THEM
def _subject_PCA(location,n_comp,save_location):
    list_of_nii = _get_list_of_ext(location, ".nii")
    length = len(list_of_nii)
    for i in range(length):
        subloc = list_of_nii[i]
        pat_img = nib.load(subloc)
        voxtime_dat = _obtain_vt_data(pat_img).T
        print("Performing PCA using " + str(n_comp) + " components for subject " + str(subloc[4:7]) + "...")
        PCA_red = _do_PCA_v2(voxtime_dat, n_comp)
        if not os.path.exists(save_location):
            os.makedirs(save_location)
        np.save(os.path.join(save_location, 'PCA_' + str(subloc[4:7] + '.npy')), PCA_red)
        print("PCA reduction done for subject " + str(subloc[4:7]) + ".")
        del pat_img, voxtime_dat, PCA_red
    print('PCA reduction -- DONE')

loc = 'E:/LRCBH/Results/1.PCA/v1'
pca_tcat = _temporal_concat(loc, vox)

import h5py
import hdf5storage

hdf5storage.write(pca_tcat,pca_result,'pca_tcat_uw.mat',matlab_compatbile='True')