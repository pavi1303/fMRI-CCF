clear all

pca_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\1.PCA\Useful';
qa_savedir ='E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\2.DR\QA';
cd('E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\1.ICA');
load('gica_100_result.mat');
dirpath = dir(dr_savedir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];

dirloc = dir(pca_savedir);
subpath = {dirloc.name}';
subpath(ismember(subpath,{'.','..'})) = [];

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

clearvars -except reordered_idx subloc subpath dr_savedir pca_savedir S;
%------------------------------------------------------------------------------%
%                 QA PROCESS - SHORTLISTING PATIENTS NEEEDING FLIP
%------------------------------------------------------------------------------%
for i= 1:length(subloc)
    fprintf('Doing the QA for subject %d..\n',i);
    sub = subloc{i};
    current = strcat(dr_savedir,'\', sub);
    cd(current);
    dr = load('dualregression.mat');
    tc = (dr.ss_tc)';
    tc_reordered = tc(:,reordered_idx);
    %sm= (dr.ss_sm)';

    suboi = subpath{i};
    current = strcat(pca_savedir,'\', suboi);
    sub_data = load(current,'vt_data');
    sub_data = (sub_data.vt_data);

    for j=1:100
        zscore=10;
        q_index=[];
        while length(q_index) <100
            q_index=find(abs(S(j,:))>zscore);
            zscore=zscore*0.95;
            if zscore < 0.5
                break;
            end
        end
        voxel_index{i,j} = q_index;
        YY=[sub_data(:,q_index) tc(:,j)];
        R=corrcoef(YY);
        L=length(q_index);
        RR(i,j)=sum(R(1:L,L+1));
    end
end
RR_OI = RR(:,1:74);
[row col] = find(RR_OI<0);
for i = 1:size(idx,1)
    val(i,1) = RR_OI(idx(i,1),idx(i,2));
end
rsn = unique(col);
pat = unique(row);
pat = subloc{pat,1};
idx = horzcat(row, col, val);
idx = sortrows(idx);
pat_freq = [pat,histc(idx(:,1),pat)];
rsn_freq = [rsn,histc(idx(:,2),rsn)];

save(fullfile('E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results',sprintf('QA_results_ss_sm_v1.mat')),...
    'reordered_idx','RR','RR_OI','idx','pat','rsn','voxel_index','val','pat_freq','rsn_freq');
%------------------------------------------------------------------------------%
%                 QA PROCESS - ACTUAL FLIPPING OF SIGNS 
%------------------------------------------------------------------------------%
% Importing the mask and other details
qa_savedir ='E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\2.DR\QA';
cd('E:\LRCBH\Data\MNI_segmented');
%cd('W:\MNI_segmented')
m = load_untouch_nii('standard_binary.nii');
M = m.img;
nii_temp = m;
[x, y, z] = size(M);
mask_temp = reshape(M, [1,x*y*z]);
[~, indices] = find(mask_temp);
% FLIPPING OF SIGNS OF THE SPATIAL MAPS AND THE TIME COURSES FOR THE
% SHORTLISTED SUBJECTS
dr_savedir ='E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\2.DR\Useful';
dirpath = dir(dr_savedir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
cd('E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results');
dat = load('QA_results_S_v1.mat');
pat = dat.pat;
idx = dat.idx;
% Finding the RSN'S that need to be changed for each of the patients
for i = 1:length(pat)
    rsn_loc{i,1} = (idx(idx(:,1)==pat(i,1),2))';
end
% Reducing the patients list to just the patients that require flipping
for j = 1:length(pat)
    subloc1{j,1} = subloc{pat(j,1),1};
end
k=1;
for k=1:length(subloc1)
    rsn_idx = rsn_loc{k,1};
    sub = subloc1{k};
    current = strcat(dr_savedir,'\', sub);
    cd(current);
    dr = load('dualregression.mat');
    ss_tc = dr.ss_tc;
    ss_sm = dr.ss_sm;
    for l=1:length(rsn_idx)
        ss_tc(rsn_idx(1,l),:)= -ss_tc(rsn_idx(1,l),:);
        ss_sm(:,rsn_idx(1,l))= -ss_sm(:,rsn_idx(1,l));
    end
    % Replacing the old values with the new adjusted ones
    ss_dr_dir = strcat(qa_savedir,'\',sub);
    if ~exist(ss_dr_dir,'dir')
        mkdir(ss_dr_dir);
    end
    save(fullfile(ss_dr_dir,sprintf('dualregression.mat')),'ss_tc',"ss_sm");
    % Saving the independent components as nii images
    save_ica_nii(ss_sm',x,y,z,indices,m,'ica_',ss_dr_dir);
end
