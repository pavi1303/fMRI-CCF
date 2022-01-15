 addpath(genpath('C:\Users\PATTIAP\Desktop\COBRE_VF\Results\MNI_segmented'))
 addpath(genpath('C:\Users\PATTIAP\Desktop\Dataset\COBRE_fMRI_MNI'))
 addpath(genpath('C:\Users\PATTIAP\Downloads'))

 % Creating the nifti mask
 c1 = niftiread('c1MNI152_T1_2mm.nii');
 c1_infor = niftiinfo('c1MNI152_T1_2mm.nii')
 c2 = niftiread('c2MNI152_T1_2mm.nii');
mask = uint8(imadd(c1,c2));
mask = uint8(imbinarize(mask));
info_func = niftiinfo('MNI_012_00001.nii')
aff = info_func.Transform.T
niftiwrite(mask,'MNI_152_mask_v2.nii',c1_infor)
info_mask = niftiinfo('MNI_152_mask_v2.nii')
infor = niftiinfo('MNI152_T1_2mm.nii')

niftiwrite(mask1,'MNI_brian_mask1.nii')
img = niftiread('MNI-008.nii');
mask = single(niftiread('MNI_152_mask_v2.nii'));
inform = niftiinfo('MNI-008.nii');
mask = double(mask);
cd('C:\Users\PATTIAP\Desktop\COBRE_VF\Results')

img_masked = bsxfun(@times,img,mask);
niftiwrite(img_masked,'MNI-008_masked.nii',inform);

niftiwrite(img_masked,'MNI_masked.nii',inform);

method=11;
a1=1;
var_normal=1;
eps=1E-6;
A0=[];
shift=[];
Sigma=0;
determine_flip=1;
npca=30;

[S,W,White,E,eigval,convergence,A,B,A_reduced,X_reduced,Sigma_reduced]=...
    ica_DC_improved_v1(pca_tcat_1,Sigma,method,eps,npca,A0,a1,var_normal,shift,determine_flip);

(X,Sigma,method,eps,npca,A0,a1,var_normal,shift)
sliceViewer(mask);
pca_tcat = vertcat(pca_tcat1,pca_tcat2,pca_tcat3,pca_tcat4);
metric_z = zscore(Metric);
mm = mean(Metric);
mm_sd = std(Metric)
gICA_30_trial = X_reduced;

info_mask= niftiinfo('MNI_masked1.nii');
info_mask.Transform.T = aff
niftiwrite(mask_img,'MNI_masked1.nii',info_mask);
info_func = niftiinfo('MNI_012_00001.nii')
X = niftired('')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files =dir(fullfile('E:\LRCBH\Results\1.PCA\Wrong_results\1.PCA\Final','*.mat'))
cd('E:\LRCBH\Results\1.PCA\Wrong_results\1.PCA\Final')
list = struct2cell(files)
list_mat = list{1, :}
v=[];
for i = 1: length(files.name)
    v1 = load(['pca_tcat_part' num2str(i) '.mat' ])
    v2 = vertcat(v,v1);
end
%%%%%%%------Doing the ICA-----------------------------%%%%%%%%%%%%
cd('C:\Users\PATTIAP\Desktop\COBRE_VF\Results\1.PCA\Final')
pca_tcat_1 = vertcat(pca_tcat1,pca_tcat2,pca_tcat3,pca_tcat4,pca_tcat5);
pca_tcat_2 = vertcat(pca_tcat6,pca_tcat7,pca_tcat8,pca_tcat9,pca_tcat10,pca_tcat11);
pca_tcat = vertcat(pca_tcat_1, pca_tcat_2);

gICA_30_vt  = double(X_reduced);
ICA_4D = reshape(gICA_30,[91,109,91]);