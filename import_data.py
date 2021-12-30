import nibabel as nib
import os
import numpy as np
from sklearn.decomposition import FastICA
# The current path of the data
# I'm not able to automate this as of now because I'm not able to get to iterate over the directories
path = 'C:/Users/PATTIAP/Desktop/Dataset/MNI dataset/sub_id_48'
if os.path.exists("C:/Users/PATTIAP/Desktop/Dataset/MNI dataset/sub_id_48"):
    os.chdir("C:/Users/PATTIAP/Desktop/Dataset/MNI dataset/sub_id_48")
else:
    print("Current working directory doesn't exist")
# Appending the fMRI time series and removing the first five volumes
newlist = []
for files in os.listdir():
    if files.endswith(".nii"):
        newlist.append(files)
    else:
        continue
for i in range(len(newlist) - 1, -1, -1):
    if newlist[i].startswith(('s', 'r')):
        del (newlist[i])
# Sorting the list of nii files in the ascending order of their names
newlist.sort()
# for i in range(5):
#   newlist.pop(0)
img_cat = nib.concat_images(newlist)
save_path = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
nib.save(img_cat, os.path.join(save_path, 'MNI-048.nii'))
del newlist, img_cat

# Loading all the 4d nii files
os.chdir('C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI')
nii_list = []
for file in os.listdir():
    if file.endswith(".nii"):
        nii_list.append(file)
# Sorting of the files names
nii_list.sort()

# Temporal concatenation of the data
length = len(nii_list)
init_time=902629
init_sources=850
all_data = np.empty([init_sources, init_time])
all_data[:] = np.nan
cnt = 0
n_voxels = []
n_trs = []
tempcat_mat = np.zeros((np.multiply(n_trs,length),902629))
tempcat_mat = np.empty([850,902629])
i=1
for i in range(length):
    pat_img = nib.load(nii_list[i])
    n_voxels = np.prod(pat_img.shape[:-1])
    n_trs = pat_img.shape[-1]
    data = pat_img.get_fdata()
    voxtime_mat = (data.reshape((n_voxels, n_trs))).transpose()
    # tempcat_mat[i] = np.std(voxtime_mat, axis=0)
    #tempcat_mat.append([voxtime_mat])
    tempcat_mat = np.append(tempcat_mat,voxtime_mat,axis=0)

tempcat_mat1 = tempcat_mat
ICA_mat = tempcat_mat1[850:5100:1]
all_data[0:n_voxels, cnt:(cnt + voxtime_mat.shape[1])] = voxtime_mat
    cnt += data.shape[1]

ICA_mat = np.concatenate((tempcat_mat[0,:],tempcat_mat[1,:]), axis=0)

# Preprocessing for ICA analysis - This is working perfectly fine
# 1.Centering of the data
def center_mean(x):
    mean = np.mean(x, axis=0, keepdims=True)
    centered = x - mean
    return centered, mean


# 2. Whitening operation based on covariance matrix
# Caclulation of the covariance matrix based on Eigen Value Decomposition
def cov(x):
    mean = np.mean(x, axis=0, keepdims=True)
    n = np.shape(x)[1] - 1
    m = x - mean
    return (m.dot(m.T)) / n


def whiten(x):
    # Calculate the covariance matrix
    coVarM = cov(x)
    # Single value decomposition
    U, S, V = np.linalg.svd(coVarM)
    # Calculate diagonal matrix of eigenvalues
    d = np.diag(1.0 / np.sqrt(S))
    # Calculate whitening matrix
    whiteM = np.dot(U, np.dot(d, U.T))
    # Project onto whitening matrix
    Xw = np.dot(whiteM, x)
    return Xw, whiteM


# Preprocessing of the signals
# Cener the signals
Xc, meanX = center_mean(ICA_mat)
# Whiten the signals
Xw, whitenM = whiten(Xc)
# Checking if the preprocessing happened right
Id_mat = np.identity(Xw.shape[0])
Cov_xW = np.round(cov(Xw))

if (abs(Cov_xW) == Id_mat).all:
    pass
else:
    print("The ICA preprocessing wasn't done right")
print(np.round(cov(Xw)))

arr1 = voxtime_mat
arr2 = voxtime_mat
# Concatenating operation
# axis = 0 implies that it is being done row-wise
arr3 = (np.concatenate((arr1, arr2), axis=1)).transpose()
ICA_mat = arr3.transpose()
del ICA_mat1

ICA_mat = voxtime_mat
del voxtime_mat

fastica = FastICA(n_components=3, whiten=False)
S_ = fastica.fit_transform(ICA_mat)



# Trial trial trial
sub_dir = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
fid = open(sub_dir,'r')