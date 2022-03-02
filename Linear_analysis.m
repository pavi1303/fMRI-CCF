clc
clear all;

fcn_dir='E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\3.FCN';
fluency_dir = 'E:\LRCBH\Projects\COBRE\Results\Documents\Excel';

%------------------------------------------------------------------------------%
%      LINEAR REGRESSION - CORRELATION WITH FLUENCY  SCORE
%------------------------------------------------------------------------------%

% GENERATION OF THE DESIGN MATRIX
cd('E:\LRCBH\Projects\COBRE\Results\Documents\Excel');
fluency_ratio = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 3 107 3]);
grp = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 4 107 4]);
age = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 6 107 6]);
ed = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 7 107 7]);
interaction = readmatrix('Cobre_fluency_study_v2.xlsx', 'Sheet','regression','Range',[2 5 107 5]);
%suvr = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 8 103 8]);
%pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 9 103 9]);
%sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 10 103 10]);
suvr_dis = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 11 107 11]);

%regressor = horzcat(fluency_ratio, grp);
regressor = horzcat(fluency_ratio, grp, interaction);
covariates = horzcat(age, ed, suvr_dis);
X = horzcat(regressor, covariates);

% GENERATION OF THE RESPONSE VARIABLE
% Getting the indices for the different RSN groups
visual_idx = [1,20,29,45,65,68,70,75,90,95];
auditory_idx = [22,36,40,48];
language_idx = [5,6,14,46,59,61,77,80,86,97,100,30,34,63];
mem_cog_idx = [16,38,50,54,56,57,60,64,82,83,89,39,41,72,94,73];
subcortical_idx = [3,9,11,13,15,42,69,74,76];
cerebellar_idx =[8,18,21,24,27,37,53,58,87,91,98];
motor_idx =[7,23,43,81,88,92];
sensory_idx=[10,26,49,19];
noise_idx=[2,4,12,17,25,28,31,32,35,44,47,51,52,55,62,66,67,71,78,79,84,85,93,96,99,33];
reordered_idx = horzcat(visual_idx,auditory_idx,language_idx,...
    mem_cog_idx,subcortical_idx,cerebellar_idx,motor_idx,sensory_idx,noise_idx);

% Reordering the time courses to generate subject-wise functional
% connectiviy matrices
dr_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\2.DR\Useful';
fcn_savedir='E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\4.FCN_reordered';
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
    % Making sure the derived functional connectivity matrix reflects the
    % required order of the different rsn's
    tc_reordered = tc(:,reordered_idx);
    [pval_reordered, fcn_reordered] = corrcoef(tc_reordered);
    save(fullfile(fcn_savedir, sprintf('FCN_reordered_%s.mat',sub)),'fcn_reordered','pval_reordered');
end

% Creating the Y - matrix
dirloc = dir(fcn_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
% Generating the patient * fcn value matrix
for j = 1:length(subloc)
    sub = subloc{j};
    current = strcat(fcn_savedir,'\',sub);
    sub_fcn = load(current,'fcn_reordered');
    sub_fcn = sub_fcn.fcn_reordered;
    sub_fcn(:,75:end) = [];
    sub_fcn(75:end,:) = [];
    sub_fcn1 = sub_fcn.';
    m  = tril(true(size(sub_fcn1)),-1);
    fcn_mat(j,:) = sub_fcn1(m).';
end
% Needs to be changed so that the matrix is patients * groups and not
% patients * functional connections
Y = fcn_mat;




