
import nibabel as nib
import os
import nilearn as nil

#The current path of the data
path = 'Z:\COBRE_SCANS/013/105780795_012_cmrr_mbep2d_bold_AP_MB8_2mm_ISO_20170908/intermediate'
if os.path.exists("Z:/COBRE_SCANS/013/105780795_012_cmrr_mbep2d_bold_AP_MB8_2mm_ISO_20170908/intermediate"):
    os.chdir("Z:/COBRE_SCANS/013/105780795_012_cmrr_mbep2d_bold_AP_MB8_2mm_ISO_20170908/intermediate")
else:
    print("Current working directory doesn't exist")
newlist=[]
for files in os.listdir():
    if files.endswith(".nii"):
        newlist.append(files)
    else:
        continue
for i in range(len(newlist)-1,-1,-1):
    if newlist[i].startswith(('s','r')):
        del(newlist[i])
newlist.sort()
for i in range(5):
    newlist.pop(0)
img_cat = nil.image.concat_imgs(newlist)
nii4D = nib.concat_images(newlist)

newlist1=[];
dir_path = 'Z:/COBRE_SCANS/008'
subdir_inc = 'intermediate'
for root, dirs, files in os.walk('Z:/COBRE_SCANS/008/'):
    if dir == 'intermediate':
        print(dir)



    for file in filelist:
        if file.endswith(".nii"):
                newlist1.append(file)
    for dir1 in subdir_inc:
        print(dir1)
os.getcwd()

del img_cat