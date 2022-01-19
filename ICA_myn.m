addpath(genpath('C:\Users\PATTIAP\Desktop\COBRE_VF\Results\MNI_segmented'))
addpath(genpath('C:\Users\PATTIAP\Documents\GitHub\fMRI-multimodal\Multivariate decomposition\ICA'))
addpath(genpath('C:\Users\PATTIAP\Desktop\COBRE_VF\Results'))
addpath(genpath('E:\LRCBH\COBRE-MNI\Individual_data\046'))

% Creating the nifti mask
c1 = niftiread('c1MNI152_T1_2mm.nii');
c1_infor = niftiinfo('c1MNI152_T1_2mm.nii')
c2 = niftiread('c2MNI152_T1_2mm.nii');
mask = uint8(imadd(c1,c2));
mask = int16(imbinarize(mask));
mask_header = niftiinfo('MNI152_T1_2mm.nii')
niftiwrite(mask,'MNI_152_mask_new.nii',mask_header);

%%% Vertical concatenation of the data
list = [8,13,16,21,23,24,27,28,29,37,38,42,43,45,46,47,48,49,50,51,52,54,55,62,66,68,69,70,72,...
    73,74,75,78,84,85,86,89,92,93,96,98,101,105,110,116,117,127,129,131,132,134,135,140,143,...
    147,150,151,153,155,157,158,159,167,169,170,174,175,176,177,178,185,187,190,191,192,...
    194,196,199,200,201,203,206,207,209,210,211,212,214,215,216,217,218,219,220,221,...
    222,223,225,228,229,230,231,233,234,235,236,237,242];
no = length(list);
location = 'E:\LRCBH\Results\1.PCA\Updated\With_standardization';
for i=1:no
    pat = list(i);
    file = fullfile(location, sprintf('PCA_%d.mat', pat));
    data = load(file);
    data = struct2cell(data);
    pat_data{1,i} = data;
end
pca_tcat = cellfun(@vertcat, pat_data);
tcat_data = vertcat(pca_tcat{:});
data = tcat_data(21600,:);
%%%%%%%------Doing the ICA-----------------------------%%%%%%%%%%%%
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
    ica_DC_improved(tcat_data,Sigma,method,eps,npca,A0,a1,var_normal,shift,determine_flip);

gICA_30_pca  = S;
ICA_4D = reshape(gICA_30,[91,109,91]);

res = isfinite(PCA_red);
res_loc = find(res==0)