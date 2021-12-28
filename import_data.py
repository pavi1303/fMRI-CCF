import nibabel as nib
import os
import numpy as np

#The current path of the data
#I'm not able to automate this as of now because I'm not able to get to iterate over the directories
path = 'C:/Users/PATTIAP/Desktop/Dataset/MNI dataset/48'
if os.path.exists("C:/Users/PATTIAP/Desktop/Dataset/MNI dataset/48"):
    os.chdir("C:/Users/PATTIAP/Desktop/Dataset/MNI dataset/48")
else:
    print("Current working directory doesn't exist")
#Appending the fMRI time series and removing the first five volumes
newlist=[]
for files in os.listdir():
    if files.endswith(".nii"):
        newlist.append(files)
    else:
        continue
for i in range(len(newlist)-1,-1,-1):
    if newlist[i].startswith(('s','r')):
        del(newlist[i])
#Sorting the list of nii files in the ascending order of their names
newlist.sort()
#for i in range(5):
 #   newlist.pop(0)
img_cat = nib.concat_images(newlist)
save_path = 'C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI'
nib.save(img_cat,os.path.join(save_path,'MNI-048.nii'))
del newlist, img_cat

#Loading all the 4d nii files
os.chdir('C:/Users/PATTIAP/Desktop/Dataset/COBRE_fMRI_MNI')
nii_list=[]
for file in os.listdir():
    if file.endswith(".nii"):
        nii_list.append(file)
#Sorting of the files names
nii_list.sort()

#Temporal concatenation of the data
length = len(nii_list)
n_voxels = []
n_trs = []
#tempcat_mat = np.zeros((length,850))
tempcat_mat =[]
for i in range(length):
    pat_img = nib.load(nii_list[i])
    n_voxels = np.prod(pat_img.shape[:-1])
    n_trs = pat_img.shape[-1]
    data = pat_img.get_fdata()
    #n_voxels.append()
    #n_trs.append()
    voxtime_mat = data.reshape((n_voxels, n_trs))
    #tempcat_mat[i] = np.std(voxtime_mat, axis=0)
    tempcat_mat.append([voxtime_mat])

#ICA_mat = tempcat_mat.transpose()

#Preprocesing for ICA analysis
#1.Centering of the data
def center_mean(x):
    mean = np.mean(x, axis=1, keepdims=True)
    centered =  x - mean
    return centered, mean
#2. Whitening operation based on covariance matrix
#Caclulation of the covariance matrix based on Eigen Value Decomposition
def cov(x):
    mean = np.mean(x, axis=1, keepdims=True)
    n = np.shape(x)[1]-1
    m = x - mean
    return (m.dot(m.T))/n

def whiten(x):
    # Calculate the covariance matrix
    coVarM = cov(X)
    # Single value decomposition
    U, S, V = np.linalg.svd(coVarM)
    # Calculate diagonal matrix of eigenvalues
    d = np.diag(1.0 / np.sqrt(S))
    # Calculate whitening matrix
    whiteM = np.dot(U, np.dot(d, U.T))
    # Project onto whitening matrix
    Xw = np.dot(whiteM, X)
    return Xw, whiteM

#Preprocessing of the signals
#Cener the signals
Xc, meanX = center_mean(arr3)
#Whiten the signals
Xw, whitenM = whiten(Xc)









arr1 = voxtime_mat
arr2 = voxtime_mat
# Concatenating operation
# axis = 0 implies that it is being done row-wise
arr3 = (np.concatenate((arr1, arr2), axis=1)).transpose()
ICA_mat = arr3.transpose()
del ICA_mat