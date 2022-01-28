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


% Generating the FCN matrix - not needed now 
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


% Generating the design matrix
fluency_ratio = readmatrix('Cobre_fluency_study.xlsx','Sheet','regression','Range',[2 3 89 3]);
% pf = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 8 89 8]);
% sf = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 9 89 9]);
age = readmatrix('Cobre_fluency_study.xlsx','Sheet','regression','Range',[2 4 89 4]);
ed = readmatrix('Cobre_fluency_study.xlsx','Sheet','regression','Range',[2 5 89 5]);
suvr = readmatrix('Cobre_fluency_study.xlsx','Sheet','regression','Range',[2 6 89 6]);
grp = readmatrix('Cobre_fluency_study.xlsx','Sheet','regression','Range',[2 7 89 7]);
regressor = horzcat(fluency_ratio, grp);
covariates = horzcat(age, ed, suvr);
interaction = fluency_ratio.*grp;
 
% Fitting the multiple linear regression model
% With the interaction term
[regress_nointer.X, regress_nointer.Y, regress_nointer.fit_model, regress_nointer.betas, regress_nointer.pvalue, ...
    regress_nointer.tstatistic, regress_nointer.alpha_level, regress_nointer.pval_interactionvar, regress_nointer.meanRsquared_original, regress_nointer.meanRsquared_adjusted, regress_nointer.sig_voxel] = ...
    regress_model('E:\LRCBH\Results\Matlab\3.DR\Unbiased',regressor,interaction, covariates,'dualregression',14,'interaction','E:\LRCBH\Results\Matlab\v2\5.Association');
% Without the interaction term
[grp_res.Yfitted,grp_res.Yfitted_sig,grp_res.beta, ...
    grp_res.pvalue, grp_res.pval_sig,grp_res.tstatistic,grp_res.alpha,grp_res.sig_asso] = confound_fitlm('E:\LRCBH\Results\Matlab\v2\4.FCN\All',X,'fcn', 19,'all',asso_savedir);

% Run and make sure it runs; slightly modify the code so that you save the
% R2 coefficient and the mean Rsquared value as well.

Y = reg_nointer.Y;
lr_model = reg_nointer.fit_model;
for k = 1:size(Y,2)
    Rsquared_orig(k,:) = fit_model{k,1}.Rsquared.Ordinary;
    Rsquared_adjust(k,:) = fit_model{k,1}.Rsquared.Adjusted;
end
Rsquared_orig_mean = mean(Rsquared_orig);
Rsquared_adjust_mean = mean(Rsquared_adjust);
