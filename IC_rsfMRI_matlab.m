%Getting the required subdir
clc
clear
tic
% Getting the list of all the directories
rootdir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data';
pca_savedir = 'E:\LRCBH\Results\Matlab\1.PCA';
ica_savedir = 'E:\LRCBH\Results\Matlab\2.ICA';
dirpath = dir(rootdir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
%Loading the mask file
cd('E:\LRCBH\MNI_segmented');
%cd('W:\MNI_segmented')
m = load_untouch_nii('standard_binary.nii');
M = m.img;
nii_temp = m;
[x, y, z] = size(M);
mask_temp = reshape(M, [1,x*y*z]);
[~, indices] = find(mask_temp);
trs = 850;
comp = 200;
%Iterating through each of the subjects
for i=1:length(subloc)
    suboi = subloc{i};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    files = dir('*.nii');
    temp = zeros(trs, 228453);
    fprintf('Generating voxeltime data for subject %s...\n',suboi);
    [temp,index] = generate_vt(files,M,x,y,z);
    vt_data = double(zscore(temp(16:end, :)));
    clear temp;
    fprintf('Performing PCA reduction using %d components for subject %s...\n',comp,suboi);
    %Performing subject wise PCA reduction
    PCA_red = double(do_PCA(vt_data,comp));
    cd(pca_savedir);
    save(fullfile(pca_savedir, sprintf('PCA_%s.mat',suboi)),'PCA_red','vt_data','index');
end
fprintf('PCA reduction done.\n');
% Doing temporal concatenation of the subjects
fprintf('Performing temporal concatenation of subjects...\n');
tcat_data = temporal_concat(subloc,pca_savedir);
fprintf('Temporal concatenation done.\n')
clearvars -except tcat_data
% Performing spatial ICA based on hyperbolic tangent
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
save(fullfile(ica_savedir,sprintf('gica_%d_result.mat',npca)),'A','S','W','White','-v7.3');
toc

%Saving the individual ica components as nii files
for j = 1:size(S,1)
    temp = zeros(1, x*y*z);
    ica_comp = S(j,:);
    temp(1,indices) = ica_comp;
    ica_3d = reshape(temp,[x,y,z]);
    nii_temp = m;
    nii_temp.img = [];
    nii_temp.img = single(ica_3d);
    Sname = ['ica_',sprintf(['%03d.nii'],j)];
    save_untouch_nii(nii_temp,[ica_savedir,'\',Sname]);
end
% Doing dual regression
cd(ica_savedir);
ica_result=load('gica_30_result.mat');
dr_savedir = 'E:\LRCBH\Results\Matlab\3.DR';
pca_savedir = 'E:\LRCBH\Results\Matlab\1.PCA\Trial';
S=ica_result.S;
noise_idx = [2,4,11,13,14,17,18,19,28,23];
%Using only the useful components
S(noise_idx,:)=[];
%Getting the directory having the subject wise voxel time data
dirloc = dir(pca_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(pca_savedir,'\', suboi);
    sub_data = load(current,'vt_data');
    sub_data = (sub_data.vt_data)';
    %Performing spatial regression
    ss_tc = (pinv(S))'*sub_data;
    %Performing temporal regression
    ss_sm = sub_data*pinv(ss_tc);
    ss_dr_dir = strcat(dr_savedir,'\',suboi(5:7));
    if ~exist(ss_dr_dir,'dir')
        mkdir(ss_dr_dir);
    end
    save(fullfile(ss_dr_dir,sprintf('dualregression.mat')),'ss_tc',"ss_sm");
    % Saving the independent components as nii images
    save_ica_nii(ss_sm',x,y,z,indices,m,ss_dr_dir);
end



    i=1;

    x = (pinv(ss_tc));
    x=x';

%List of inputs to the dual regression function
% The group ICA map
% the noise components indices
% dual regression save dir
% pca directory (voxel time data files)
heatmap(R);
