 addpath(genpath('C:\Users\PATTIAP\Desktop\COBRE_VF\Results\MNI_segmented'))
 addpath(genpath('C:\Users\PATTIAP\Desktop\Dataset\COBRE_fMRI_MNI'))
 addpath(genpath('C:\Users\PATTIAP\Desktop\COBRE_VF\Results'))
 addpath(genpath('E:\LRCBH\COBRE-MNI\Individual_data\046'))

 % Creating the nifti mask
 c1 = niftiread('c1MNI152_T1_2mm.nii');
 c1_infor = niftiinfo('c1MNI152_T1_2mm.nii')
 c2 = niftiread('c2MNI152_T1_2mm.nii');
mask = uint8(imadd(c1,c2));
mask = int16(imbinarize(mask));
M1 = max(mask, [], 'all');
M2 = min(mask, [], 'all');
info_func = niftiinfo('MNI_012_00001.nii')
mask_header = niftiinfo('MNI152_T1_2mm.nii')
aff = info_func.Transform.T
niftiwrite(mask,'MNI_152_mask_new.nii',mask_header);
mask_new = niftiread('MNI_152_mask_new.nii');
M1 = max(mask_new, [], 'all');
M2 = min(mask_new, [], 'all');
info_mask = niftiinfo('MNI_152_mask_v2.nii')
infor = niftiinfo('MNI152_T1_2mm.nii')
%%%% Creating the NIFTI mask
img = niftiread('s8_MNI_013_00001.nii');
niftiwrite(mask1,'MNI_brian_mask1.nii')
img = niftiread('MNI-008.nii');
mask = single(mask_new);
mask = single(niftiread('MNI_152_mask_v2.nii'));
inform = niftiinfo('s8_MNI_013_00001.nii');
mask = double(mask);
cd('C:\Users\PATTIAP\Desktop\COBRE_VF\Results')

img_masked = bsxfun(@times,img,mask);
niftiwrite(img_masked,'matlab_masked.nii',inform);

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
%%% Vertical concatenation of the data
list = [8,13,16,21,23,24,27,28,29,37,38,42,43,45,46,47,48,49,50,51,52,54,55,62,66,68,69,70,72,...
    73,74,75,78,84,85,86,89,92,93,96,98,101,105,110,116,117,127,129,131,132,134,135,140,143,...
    147,150,151,153,155,157,158,159,167,169,170,174,175,176,177,178,185,187,190,191,192,...
    194,196,199,200,201,203,206,207,209,210,211,212,214,215,216,217,218,219,220,221,...
    222,223,225,228,229,230,231,233,234,235,236,237,242];
no = length(list);
location = 'E:\LRCBH\Results\1.PCA\Updated\No_standardization';
for i=1:no
    pat = list(i);
    file = fullfile(location, sprintf('PCA_%d.mat', pat));
    data = load(file);
    data = struct2cell(data);
    pat_data{1,i} = data;
end
pca_tcat = cellfun(@vertcat, pat_data);
tcat_data = vertcat(pca_tcat{:});
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
S1 = S;

gICA_30_vt_v2  = S1;
ICA_4D = reshape(gICA_30,[91,109,91]);

res = isfinite(PCA_red);
res_loc = find(res==0)