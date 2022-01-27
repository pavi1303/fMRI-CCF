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
asso_savedir = 'E:\LRCBH\Results\Matlab\v2\5.Association';
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

% Controlling for confounds and finding significant associations
% Generating the design matrix
fluency_ratio = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 3 89 3]);
pf = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 8 89 8]);
sf = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 9 89 9]);
age = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 4 89 4]);
ed = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 5 89 5]);
suvr = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 6 89 6]);
grp = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 7 89 7]);
X = horzcat(fluency_ratio, age, ed, suvr, grp);
X1 = horzcat(fluency_ratio, age, ed, suvr);
X_grp1 = X1(1:44,:);
X_grp2 = X1(45:end,:);

% Fitting the linear regression model 
[grp1_res.Yfitted,grp1_res.Yfitted_sig,grp1_res.beta, ...
    grp1_res.pvalue, grp1_res.pval_sig, grp1_res.tstatistic,grp1_res.alpha,grp1_res.sig_asso] = confound_fitlm('E:\LRCBH\Results\Matlab\v2\4.FCN\Grp1',X_grp1,'fcn', 19,'grp1',asso_savedir);
[grp2_res.Yfitted,grp2_res.Yfitted_sig,grp2_res.beta,...
    grp2_res.pvalue, grp2_res.pval_sig,grp2_res.tstatistic,grp2_res.alpha,grp2_res.sig_asso] = confound_fitlm('E:\LRCBH\Results\Matlab\v2\4.FCN\Grp2',X_grp2,'fcn', 19,'grp2',asso_savedir);
[grp_res.Yfitted,grp_res.Yfitted_sig,grp_res.beta, ...
    grp_res.pvalue, grp_res.pval_sig,grp_res.tstatistic,grp_res.alpha,grp_res.sig_asso] = confound_fitlm('E:\LRCBH\Results\Matlab\v2\4.FCN\All',X,'fcn', 19,'all',asso_savedir);

% Generating the scatter plot (VF vs FC)
VF = repmat(fluency_ratio,[1 16]);
Y = grp_res.Yfitted_sig;
g = horzcat(repmat(1,[1,44]),repmat(0,[1,44]))';
h = gscatter(VF,Y,g);
gscatter(fluency_ratio,Y(:,1),g)

%
[f1,xi1] = ksdensity(pf);
[f2,xi2] = ksdensity(sf);
[f3,xi3] = ksdensity(fluency_ratio);

figure;
plot(xi1,f1)
lgd = legend('Phonemic fluency');
title('pd estimate - Phonemic fluency');
figure;
plot(xi2,f2)
lgd = legend('Semantic fluency');
title('pd estimate - Semantic fluency');
figure;
plot(xi3,f3)
lgd = legend('Fluency ratio');
title('pd estimate - Fluency ratio');

pval = 0.05;
pval_adjusted = 0.05/171;