%Getting the required subdir
clc
clear
tic
% Getting the list of all the directories
rootdir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data';
pca_savedir = 'E:\LRCBH\Results\Matlab\1.PCA';
ica_savedir = 'E:\LRCBH\Results\Matlab\2.ICA';
dr_savedir = 'E:\LRCBH\Results\Matlab\3.DR';
fcn_savedir='E:\LRCBH\Results\Matlab\4.FCN';
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
ica_dat = load('gICA_30_result.mat','S');
S = ica_dat.S;
save_ica_nii(S,x,y,z,indices,m,ica_savedir)

% Doing dual regression
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
    fprintf('Performing dual regression for subject %s...\n',suboi(5:7));
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
fprintf('Dual regression done.\n')


% Generating the FCN matrix
dirpath = dir(dr_savedir);
subdir = [dirpath(:).isdir];
subpath = {dirpath(subdir).name}';
subpath(ismember(subpath,{'.','..'})) = [];
for i=1:length(subpath)
    sub = subpath{i};
    current = strcat(dr_savedir,'\', sub);
    cd(current);
    dr = load('dualregression.mat','ss_tc');
    tc = (dr.ss_tc)';
    fcn = corrcoef(tc);
    save(fullfile(fcn_savedir, sprintf('FCN_%s.mat',sub)),'fcn');
end

% Controlling for confounds
% Getting the directory having the subject wise FNC matrix data
fcn_savedir = 'E:\LRCBH\Results\Matlab\v2\4.FCN\Grp1';
%fcn_Savedir = 'E:\LRCBH\Results\Matlab\v2\4.FCN\Grp2';

% Generating the design matrix
fluency = readmatrix('EXCEL_Language_Control_Study_2021-12-03-023219028 (Autosaved).xlsx','Sheet','regression','Range',[2 3 89 3]);
age = readmatrix('EXCEL_Language_Control_Study_2021-12-03-023219028 (Autosaved).xlsx','Sheet','regression','Range',[2 4 89 4]);
ed = readmatrix('EXCEL_Language_Control_Study_2021-12-03-023219028 (Autosaved).xlsx','Sheet','regression','Range',[2 5 89 5]);
suvr = readmatrix('EXCEL_Language_Control_Study_2021-12-03-023219028 (Autosaved).xlsx','Sheet','regression','Range',[2 6 89 6]);
%grp = readmatrix('EXCEL_Language_Control_Study_2021-12-03-023219028 (Autosaved).xlsx','Sheet','regression','Range',[2 7 89 7]);
X = horzcat(fluency, age, ed, suvr);
x = ones(size(X,1),1);
X = [X x];
X_grp1 = X(1:44,:);

dirloc = dir(fcn_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(fcn_savedir,'\', suboi);
    sub_data = load(current,'fcn');
    sub_data = (sub_data.fcn)';
    corr_mat = tril(sub_data,-1);
    corr_vec = nonzeros(corr_mat);
    fnc_val{i,1} = corr_vec';
end
% Forming the Y data
Y = vertcat(fnc_val{:});
% Performing regression to obtain coefficients
for k=1:size(Y,2)
    lr_model{k,1} = fitlm(X,Y(:,k));
end

for k = 1:size(Y,2)
     [b(:,k),~,~,~,stats(:,k)] = regress(Y(:,k),X_grp1);
 end
p-val
[b1,~,~,~,stats1] = regress(Y(:,1),X_grp1);
stats(3,:) = stats(3,:)/size(Y,2);
p = 0.05/size(Y,2);
idx = find(stats(3,:)<p);

%
[row,col] = find(corr_mat);
idx1 = horzcat(row,col);
sig_idx = idx1(idx,:);
% Inputs for the confound function
% Directory of the functional connectivity matrices
% The design matrix
[coeff,pval_corr,stats_corrected,n_compar,sig_loc]
S_grp1 =struct;
S_grp2 = struct;
[S_grp1.coeff, S_grp1.pval, S_grp1.stats, S_grp1.comparisons, S_grp1.sig_asso] = confound_sig(fcn_savedir,X_grp1,'fcn',corr_mat);
%  Find the significant associations with and without the inclusion of
%  cerebellar networks

