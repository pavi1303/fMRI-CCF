clc;
clear;

% Loading the X variables
cd('F:\LRCBH\Projects\COBRE\Results\Documents\Excel');
grp = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 4 90 4]);
idx = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 17 90 17]);
pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 9 90 9]);
sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 11 90 11]);
age = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 6 90 6]);
ed = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 7 90 7]);
gen = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 16 90 16]);

X = [zscore(sf), zscore(pf), zscore(age), zscore(ed), gen];

% Loading the masks and the ROI's
cd('F:\LRCBH\Data\MNI_segmented');
gm_mask = double(load_untouch_nii('GM_binary.nii').img);
cd('F:\LRCBH\Projects\COBRE\Results\Laterality_index\ROI\IFG');
lifg_roi = double(load_untouch_nii('LIFG_binary.nii').img);
rifg_roi = double(load_untouch_nii('RIFG_binary.nii').img);
cd('F:\LRCBH\Projects\COBRE\Results\Laterality_index\ROI\STG');
lstg_roi = double(load_untouch_nii('LSTG_binary.nii').img);
rstg_roi = double(load_untouch_nii('RSTG_binary.nii').img);


timeseries_dir = 'F:\LRCBH\Data\COBRE_FLUENCY\Preprocessed timeseries';
trs_savedir = 'F:\LRCBH\Projects\COBRE\Results\Laterality_index\Timeseries';
fc_savedir = 'F:\LRCBH\Projects\COBRE\Results\Laterality_index\FC and LI';

%seed = 'IFG';
seed = 'STG';
LI_index(timeseries_dir, seed, lstg_roi, rstg_roi, lifg_roi, rifg_roi, gm_mask, trs_savedir, fc_savedir);


% Forming the laterality index vector
pat_list = dir(fullfile(fc_savedir, seed));
LI_vector = zeros(length(pat_list)-2, 1);
cd(fullfile(fc_savedir, seed));
for i = 3:length(pat_list)
    LI_vector(i-2) = load((pat_list(i).name)).LI;
end
LI_vector_new = LI_vector(idx,:);

% Vector to save for boxplot
LI_result = horzcat(LI_vector_new, grp);
save(fullfile('F:\LRCBH\Projects\COBRE\Results\Laterality_index', sprintf('LI_result_%s.mat', seed)), "LI_result");
temp = strings(length(pat_list)-2, 1);

for j = 1:length(LI_result)
    if LI_result(j,1)>0.2
        temp(j) = 'Left';
    elseif LI_result(j,1)<-0.2
        temp(j) = 'Right';
    else
        temp(j) = 'Bilateral';
    end
end

left_nc  = nnz(strcmp(temp(1:44),'Left'));
right_nc  = nnz(strcmp(temp(1:44),'Right'));
bi_nc  = nnz(strcmp(temp(1:44),'Bilateral'));

left_mci  = nnz(strcmp(temp(45:end),'Left'));
right_mci  = nnz(strcmp(temp(45:end),'Right'));
bi_mci  = nnz(strcmp(temp(45:end),'Bilateral'));

% Variance inflation factor
R0 = corrcoef(X); % correlation matrix
VIF1 =diag(inv(R0))';

[~,p,~,stats] = ttest2(LI_vector_new(find(grp==0)), LI_vector_new(find(grp==1)));

% Linear regression - Association analysis
%NC
level1_res_nc = fitlm(X(1:44,3:end), LI_result(1:44,1)).Residuals.Raw;
level2_model_nc = fitlm(X(1:44,1:2), level1_res_nc);
%MCI
level1_res_mci = fitlm(X(45:end,3:end), LI_result(45:end,1)).Residuals.Raw;
level2_model_mci = fitlm(X(45:end,1:2), level1_res_mci);


%%%
imagesc(timeseries); colormap('gray');

%%Seed based correlation analysis
timeseries_dir = 'W:\LRCBH\COBRE_rsfMRI\Useful\Timeseries';
cd('W:\LRCBH\COBRE_rsfMRI');
MNItemplate = load_untouch_nii('standard_binary.nii');
mask = double(MNItemplate.img);
sbc_savedir = 'W:\LRCBH\Laterality index\SBC maps';

cd('W:\LRCBH\Laterality index\ROI\IFG');%IFG seeds
sbc_maps(timeseries_dir, double(load_untouch_nii('LIFG_binary.nii').img), mask, MNItemplate, sbc_savedir, 'LIFG');%LIFG
cd('W:\LRCBH\Laterality index\ROI\IFG');%IFG seeds
sbc_maps(timeseries_dir, double(load_untouch_nii('RIFG_binary.nii').img), mask, MNItemplate, sbc_savedir, 'RIFG');%RIFG

cd('W:\LRCBH\Laterality index\ROI\STG');%STG seeds
sbc_maps(timeseries_dir, double(load_untouch_nii('LSTG_binary.nii').img), mask, MNItemplate, sbc_savedir, 'LSTG');%LSTG
cd('W:\LRCBH\Laterality index\ROI\STG');%STG seeds
sbc_maps(timeseries_dir, double(load_untouch_nii('RSTG_binary.nii').img), mask, MNItemplate, sbc_savedir, 'RSTG');%RSTG


%%% LI based on the spatial ICA maps %%%

% Generating thresholded ICA spatial maps
cd('F:\LRCBH\Data\MNI_segmented');
M = double(load_untouch_nii('standard_binary.nii').img);
[x, y, z] = size(M);
[~, indices] = find(reshape(M, [1,x*y*z]));
% Loading the ICA results
cd('C:\Users\PATTIAP\Dropbox\ICA\ica');
S = double(load('gICA_30_result.mat').S);
ica_thresh_savedir = 'C:\Users\PATTIAP\Dropbox\ICA\ica\thresholded';

% Creating binary thresholded maps
thr = 3;
for i = 1:size(S, 1)
    for j = 1:size(S, 2)
        if S(i, j)<thr
            S(i, j) = 0;
        elseif S(i, j)>thr
            S(i, j) = 1;
        end
    end
end
save_ica_nii(S,x,y,z,indices,m,'gICA_',ica_thresh_savedir);
% Calculating the dice similarity coefficient between the ICA thresholded
% maps and the language fROI from the Willard Atlas
% Loading the Willard Atlas
cd('F:\LRCBH\Atlas\Willard shirer atlas\fROIs_90\1.Merged');
lang_atlas = double(load_untouch_nii('Language.nii').img);
lang_atlas = lang_atlas.*M;
cd(ica_thresh_savedir);
ica_list = dir('*.nii');
dice_ica_merged = zeros(1, length(ica_list));
jaccard_ica_merged = zeros(1, length(ica_list));
for ica_comp = 1:length(ica_list)
    temp = double(load_untouch_nii(ica_list(ica_comp).name).img);
    dice_ica_merged(1, ica_comp) = dice(lang_atlas, temp);
    jaccard_ica_merged(1, ica_comp) = jaccard(lang_atlas, temp);
end
% LI calculation
cd('F:\LRCBH\Data\MNI_segmented');
MNItemplate = load_untouch_nii('standard_binary.nii');
LH_mask = double(load_untouch_nii('LH_standard_binary.nii').img);
RH_mask = double(load_untouch_nii('RH_standard_binary.nii').img);
dr_dir = 'F:\LRCBH\Projects\COBRE\Results\Laterality_index\ICA\Dual regression\Spatial maps';
cd('F:\LRCBH\Projects\COBRE\Results\Documents\Excel');
order = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 17 90 17]);
grp = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 4 90 4]);
pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 10 90 10]);
sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 11 90 11]);
age = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 6 90 6]);
ed = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 7 90 7]);

%ICA component 20
ica_idx = 20;
thr = 2;
LI_thr = 0.1;
LI_struct_20_2_01 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir);
thr = 3;
LI_struct_20_3_01 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir);
LI_thr = 0.2;
LI_struct_20_3_02 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir);
%ICA component 4
ica_idx = 4;
thr = 2;
LI_thr = 0.1;
LI_struct_4_2_01 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir);
thr = 3;
LI_struct_4_3_01 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir);
LI_thr = 0.2;
LI_struct_4_3_02 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir);
%ICA component 9
ica_idx = 9;
thr = 2;
LI_thr = 0.1;
LI_struct_9_2_01 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir);
thr = 3;
LI_struct_9_3_01 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir, S);
LI_thr = 0.2;
LI_struct_9_3_02 = LI_ICA(MNItemplate, ica_idx, thr, LH_mask, RH_mask, LI_thr, order, grp, dr_dir);

% With SF
figure;
scatter(LI_struct_9_3_01.LI_values(1:44), sf(1:44));
r_nc = corr(LI_struct_9_3_01.LI_values(1:44), sf(1:44), 'type', 'Pearson');
figure;
scatter(LI_struct_9_3_01.LI_values(45:end), sf(45:end));
r_mci = corr(LI_struct_9_3_01.LI_values(45:end), sf(45:end), 'type', 'Pearson');
% With age
r_nc = corr(LI_struct_9_3_01.LI_values(1:44), age(1:44), 'type', 'Pearson');
r_mci = corr(LI_struct_9_3_01.LI_values(45:end), age(45:end), 'type', 'Pearson');
% With education
r_nc = corr(LI_struct_9_3_01.LI_values(1:44), ed(1:44), 'type', 'Pearson');
r_mci = corr(LI_struct_9_3_01.LI_values(45:end), ed(45:end), 'type', 'Pearson');

gscatter(sf, age, grp);
r_nc = corr(zscore(sf(1:44)), zscore(age(1:44)), 'type', 'Pearson');
r_mci = corr(zscore(sf(45:end)), zscore(age(45:end)), 'type', 'Pearson');

age_rm = age(1:44);
sf_rm = sf(1:44);
mdl = fitlm(sf(1:44), age(1:44));
plotDiagnostics(mdl,'cookd');

%%%% REDUNDANT CODE
% % Getting the LH voxel indices
% cd('F:\LRCBH\Data\MNI_segmented');
% LH_standard = double(load_untouch_nii('LH_standard_binary.nii').img);
% [~, LH_indices] = find(reshape(LH_standard, [1,x*y*z]));
% % Getting the RH voxel indices
% cd('F:\LRCBH\Data\MNI_segmented');
% RH_standard = double(load_untouch_nii('RH_standard_binary.nii').img);
% [~, RH_indices] = find(reshape(RH_standard, [1,x*y*z]));
% % Obtaining the LI values - Test run
% cd(ica_thresh_savedir);
% temp = double(load_untouch_nii('gICA_009.nii').img);
% temp_trs = reshape(temp, [1,x*y*z]);
% temp_trs_LH = nnz(temp_trs(:,LH_indices_unique));
% temp_trs_RH = nnz(temp_trs(:,RH_indices_unique));
% LI_value = (temp_trs_LH-temp_trs_RH)/(temp_trs_LH+temp_trs_RH);
% %%%%%%%
% LH_indices_unique = setdiff(LH_indices, RH_indices);
% RH_indices_unique = setdiff(RH_indices, LH_indices);
% total_indices = unique(horzcat(LH_indices, RH_indices));
% % Trying for the individual language ROI's
% lang_atlas_dir = 'F:\LRCBH\Atlas\Willard shirer atlas\fROIs_90\2.Individual\Language\Original';
% cd(lang_atlas_dir);
% langatlas_list = dir('*.nii');
% cd(ica_thresh_savedir);
% ica_list = dir('*.nii');
% dice_ica_ind = zeros(length(langatlas_list), length(ica_list));
% for langatlas_comp = 1:length(langatlas_list)
%     cd(lang_atlas_dir);
%     langatlas_oi = double(load_untouch_nii(langatlas_list(langatlas_comp).name).img);
%     for ica_comp = 1:length(ica_list)
%         cd(ica_thresh_savedir);
%         ica_comp_oi = double(load_untouch_nii(ica_list(ica_comp).name).img);
%         dice_ica_ind(langatlas_comp, ica_comp) = dice(langatlas_oi, ica_comp_oi);
%     end
% end
%
%
%
%
% % DSC calculation
% template = reshape(lang_atlas, [1,x*y*z]);
% template = template(:, indices);
% img = S(20, :);
% img(img<0)=0;
% intersect = nnz(template.*img);
% tsum = nnz(template);
% isum = nnz(img);
% dice_val = (2*intersect)/(tsum*isum);


%%%%%%%%-----MASS UNIVARIATE ANALYSIS-----%%%%%%%%%%%
dr_timeseries = 'F:\LRCBH\Projects\COBRE\Results\Laterality_index\ICA\Dual regression\Timeseries';
pat_list = dir(fullfile(dr_timeseries));
noise_idx = [13, 14, 16, 18, 19, 22, 24, 25, 27];
n_ica = 30;
new_ica = setdiff([1:30], noise_idx);
fc_mat = zeros(length(pat_list), ((n_ica-length(noise_idx))*(n_ica-length(noise_idx)-1))/2);
cd('F:\LRCBH\Projects\COBRE\Results\Documents\Excel'); % Getting the values of the all the independent variables
order = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 17 90 17]);
%Variables of interest
grp = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 4 90 4]);
pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 10 90 10]);
sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 11 90 11]);
%Covariates
age = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 6 90 6]);
ed = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 7 90 7]);
gen = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 16 90 16]);

for sub = 3:length(pat_list)
    cd(dr_timeseries);
    % Loading the subject timeseries
    sub_trs = load(pat_list(sub).name).ss_tc;
    sub_trs(noise_idx, :) = [];
    fc = corrcoef(sub_trs');
    % Replace this with using the ltril instead to match the indices
    % identified using the function later
    [~, ~, fc_mat(sub-2, :)] = find(tril(fc, -1));
%     At = triu(fc, 1).';
%     m  = tril(true(size(At)), -1);
%     fc_vec  = At(m).'; % Row vector of the correlation coefficient values
%     fc_mat(sub-2, :) = fc_vec;
end
fc_mat = fc_mat(1:length(pat_list)-2,:); %Removing the unwanted entries
fc_mat = fc_mat(order, :); % Each column will be part of the design matrix
% Performing linear regression for all these connections
cov_mat = [age, ed, gen];%Design matrix
Y = sf;
% Function to return the results of linear regression in presentable format
function mlr(cov_mat, Y, fc_mat, grp, sf, pf)
if sf==1
    Y =sf;
    response_name = 'SF';
elseif pf==1
    Y =pf;
    response_name = 'PF';
end
% Type 1 - Traditional regression model Var
pred_names = {'Age', 'Education', 'Gender', 'Group', 'FC', 'Interaction'};
S = struct;
for j = 1:size(fc_mat, 2)
    X = [cov_mat, grp, fc_mat(:, j), grp.*fc_mat(:,j)];
    S.X = X;
    mdl = fitlm(X, Y, 'ResponseVar',response_name, 'PredictorVars',pred_names);
    S.adj_r2(j,1) = mdl.Rsquared.Adjusted;
    pvalue(j,:) = mdl.Coefficients.pValue;
    tstat(j,:) = mdl.Coefficients.tStat;
    S.Yfitted(:,j) = mdl.Fitted;
    S.raw_res(:,j) = mdl.Residuals.Raw;
end
S.pvalue = array2table(pvalue,'VariableNames', [{'Intercept'}, pred_names]);
S.tstat = array2table(tstat,'VariableNames', [{'Intercept'}, pred_names]);
% Getting the indices of the significant p-values if any
% For a 30*30 ica matrix
new_ica = setdiff([1:n_ica], noise_idx);
[loc(:,1), loc(:, 2)] = find(tril(rand(length(new_ica)), -1));
sig_loc = new_ica(loc(find(pvalue(:, end)<0.05), :));





% Type 2 - Two level regression model
pred_names = {'Group', 'FC', 'Interaction'};
cov_mdl = fitlm(cov_mat, Y); % Regressing out the effects of covariates
Y1 = cov_mdl.Residuals.Raw;
S1 = struct;
for j = 1:size(fc_mat, 2)
    X = [grp, fc_mat(:, j), grp.*fc_mat(:,j)];
    S1.X = X;
    mdl = fitlm(X, Y1, 'ResponseVar',response_name, 'PredictorVars',pred_names);
    S1.adj_r2(j,1) = mdl.Rsquared.Adjusted;
    S1.pvalue(j,:) = mdl.Coefficients.pValue;
    S1.tstat(j,:) = mdl.Coefficients.tStat;
    S1.Yfitted(:,j) = mdl.Fitted;
    S1.raw_res(:,j) = mdl.Residuals.Raw;
end






% 3 types of regression models implemented
% 1. Single level - have all the variables together
% 2. Two level - Regress out the covariates and then do the 2nd level with
% the remaining variables
% 3. The first two versions are grouped models - 3rd one is a separate
% model
% Do the inner function as a separate function to return all the values
% also get the connections for which p-value is less (see how to get this
% from the ICA component matrix)
cov_mdl = fitlm(cov_mat, Y); % Regressing out the effects of covariates
Y1 = cov_mdl.Residuals.Raw;
for j = 1:size(fc_mat, 2)
    X = [grp, fc_mat(:, j), grp.*fc_mat(:,j)];
    mdl = fitlm(X, Y1);
    adj_r2(j,1) = mdl.Rsquared.Adjusted;
    pvalue(j,:) = mdl.Coefficients.pValue;
    tstat(j,:) = mdl.Coefficients.tStat;
    Yfitted(:,j) = mdl.Fitted;
    raw_res(:,j) = mdl.Residuals.Raw;
end

end

% Getting the significant connections
sig_loc1 = new_ica(loc_idx);






