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
tempcat_mat = np.zeros((length,850))

for i in range(length):
    pat_img = nib.load(nii_list[i])
    n_voxels = np.prod(pat_img.shape[:-1])
    n_trs = pat_img.shape[-1]
    data = pat_img.get_fdata()
    #n_voxels.append()
    #n_trs.append()
    voxtime_mat = data.reshape((n_voxels, n_trs))
    tempcat_mat[i] = np.std(voxtime_mat, axis=0)

del pat_img, data, n_voxels, n_trs, nii_list
del tempcat_mat
del length

