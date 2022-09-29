%------------------------------------------------------------------------------%
%         DATA PREPARATION - WILLARD ATLAS BASED ANALYSIS
%------------------------------------------------------------------------------%
clc;
clear all;    
rootdir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data\Useful';
atlas_loc = 'E:\LRCBH\Atlas\Willard_with_overlap_2mm\1.Merged';
data_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\1.Data_merged';
fcn_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\3.FCN\2.Individual\1.Full_correlation';
pcn_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\3.FCN\2.Individual\2.Partial_correlation';
kcn_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\3.FCN\1.Merged\3.KendallTau';
mi_savedir_own = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\4.MI\1.Merged\MAT files\Own';
mi_savedir_inbuilt = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\4.MI\1.Merged\MAT files\Inbuilt';

% Getting the list of all the subjects
dirpath = dir(rootdir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
% Iterating throught the different functional regions
cd(atlas_loc);
files = dir('*.nii');
region_list = extractfield(files,'name')';
% GENERATION OF THE REGION * TIME SERIES DATA
for j = 1:length(subloc)
    suboi = subloc{j};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    filelist = dir('*.nii');
    timeseries_voxel = [];
    timeseries_region = [];
    fprintf('Generating data for subject %s...\n',suboi);
    for i = 1:length(region_list)
        cd(atlas_loc);
        region = load_untouch_nii(region_list{i,1});
        roi = single(region.img);
        [x, y, z] = size(roi);
        roi_temp = reshape(roi, [1,x*y*z]);
        [~, indices] = find(roi_temp);
        for k=1:length(filelist)
            cd(current);
            vol = load_untouch_nii(filelist(k).name);
            I = vol.img;
            I_R = I.*roi;
            vt = reshape(I_R, [1,x*y*z]);
            %[~,idx] = find(vt);
            vt = vt(:,indices);
            data(k,:) = vt;
            avg_data = mean(data,2);
        end
        timeseries_voxel = double(horzcat(timeseries_voxel,data));
        timeseries_region = double(horzcat(timeseries_region,avg_data));
        voxel_loc{i,1} = indices;
        clear vol I I_R vt data avg_data region roi x y z roi_temp indices;
    end
    save(fullfile(data_savedir, sprintf('ROI_%s.mat',suboi)),'timeseries_voxel','timeseries_region','voxel_loc');
    clear timeseries_region timeseries_voxel voxel_loc;
end

%------------------------------------------------------------------------------%
%      LINEAR REGRESSION - CORRELATION WITH FLUENCY  SCORE
%------------------------------------------------------------------------------%

% GENERATION OF THE DESIGN MATRIX
cd('E:\LRCBH\Projects\COBRE\Results\Documents\Excel');
fluency_ratio = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 3 107 3]);
grp = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 4 107 4]);
age = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 6 107 6]);
ed = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 7 107 7]);
%interaction = readmatrix('Cobre_fluency_study_v2.xlsx', 'Sheet','regression','Range',[2 5 107 5]);
%suvr = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 8 103 8]);
pf = zscore(readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 9 107 9]));
sf = zscore(readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 10 107 10]));
suvr_dis = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 11 107 11]);

%regressor = horzcat(fluency_ratio, grp);
%regressor = horzcat(fluency_ratio, grp, interaction);
interaction = grp.*sf;
regressor = horzcat(sf, grp, interaction);
covariates = horzcat(pf,age, ed, suvr_dis);
X = horzcat(regressor, covariates);
% Design matrices for the two groups
X_nc = X(1:51,:);
X_nc(:,2:3)=[];
X_mci = X(52:end,:);
X_mci(:,2:3)=[];

% CONNECTIVITY MATRIX (Full correlation)
dirpath = dir(data_savedir);
subpath = {dirpath.name}';
subpath(ismember(subpath,{'.','..'})) = [];
for i=1:length(subpath)
    sub = subpath{i};
    current = strcat(data_savedir,'\', sub);
    tc_data = load(current,'timeseries_region');
    tc = (tc_data.timeseries_region);
    [fcn_roi, pval] = corrcoef(tc);
    save(fullfile(fcn_savedir, sprintf('FCN_merged_%s',sub)),'fcn_roi','pval');
end
% CONNECTIVITY MATRIX (Partial correlation)
dirpath = dir(data_savedir);
subpath = {dirpath.name}';
subpath(ismember(subpath,{'.','..'})) = [];
for i=1:length(subpath)
    sub = subpath{i};
    current = strcat(data_savedir,'\', sub);
    tc_data = load(current,'timeseries_region');
    tc = (tc_data.timeseries_region);
    pcorr_roi = partialcor(tc);
    save(fullfile(pcn_savedir, sprintf('PCORR_merged_%s',sub)),'pcorr_roi');
end

% CONNECTIVITY MATRIX (Full correlation - Kendall's Tau)
dirpath = dir(data_savedir);
subpath = {dirpath.name}';
subpath(ismember(subpath,{'.','..'})) = [];
for i=1:length(subpath)
    sub = subpath{i};
    current = strcat(data_savedir,'\', sub);
    tc_data = load(current,'timeseries_region');
    tc = (tc_data.timeseries_region);
    kcorr_roi = corr(tc, 'type', 'Kendall');
    save(fullfile(kcn_savedir, sprintf('KendallTau_%s',sub)),'kcorr_roi');
end
% GENERATION OF THE RESPONSE VARIABLE (Full correlation)
dirloc = dir(fcn_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for j = 1:length(subloc)
    sub = subloc{j};
    current = strcat(fcn_savedir,'\',sub);
    sub_fcn = load(current,'fcn_roi');
    sub_fcn = abs(sub_fcn.fcn_roi);
    [Y_full(j,:), Y_full_avg(j,:)] = getfc_all(sub_fcn);
end

% GENERATION OF THE RESPONSE VARIABLE (Mutual information - Own)
dirloc = dir(mi_savedir_own);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for j = 1:length(subloc)
    sub = subloc{j};
    current = strcat(mi_savedir,'\',sub);
    sub_fcn = load(current,'MI_default');
    sub_fcn = abs(sub_fcn.MI_default);
    [Y_mi(j,:), Y_mi_avg(j,:)] = getfc_all(sub_fcn);
end
% GENERATION OF THE RESPONSE VARIABLE (Mutual information - Inbuilt)
dirloc = dir(mi_savedir_inbuilt);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for j = 1:length(subloc)
    sub = subloc{j};
    current = strcat(mi_savedir_inbuilt,'\',sub);
    sub_fcn = load(current,'MI');
    sub_fcn = abs(sub_fcn.MI);
    [Y_mi(j,:), Y_mi_avg(j,:)] = getfc_all(sub_fcn);
end
% GENERATION OF THE RESPONSE VARIABLE (Partial correlation)
dirloc = dir(pcn_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for j = 1:length(subloc)
    sub = subloc{j};
    current = strcat(pcn_savedir,'\',sub);
    sub_fcn = load(current,'pcorr_roi');
    sub_fcn = abs(sub_fcn.pcorr_roi);
    [Y_partial(j,:), Y_partial_avg(j,:)] = getfc_all(sub_fcn);
end

clearvars -except X X_mci X_nc Y_full Y_full_avg Y_partial Y_partial_avg;
clearvars -except X_nc X_mci Y_full Y_partial;
clearvars -except X X_mci X_nc Y_mi Y_mi_avg;
% Fitting the regression model
mlr_roi_model= struct;
for k=1:size(Y_mi_avg,2)
    fprintf('Fitting the linear regression model for rsn group %d...\n',k);
    mlr_roi_model.lr_model{k,1} = fitlm(X,Y_mi_avg(:,k),'RobustOpts','ols');
    mlr_roi_model.coeff(k,:) = mlr_roi_model.lr_model{k,1}.Coefficients.Estimate;
    mlr_roi_model.pval(k,:) = mlr_roi_model.lr_model{k,1}.Coefficients.pValue;
    mlr_roi_model.tstat(k,:) = mlr_roi_model.lr_model{k,1}.Coefficients.tStat;
    mlr_roi_model.Rsquared_orig(k,:) = mlr_roi_model.lr_model{k,1}.Rsquared.Ordinary;
    mlr_roi_model.Rsquared_adjust(k,:) = mlr_roi_model.lr_model{k,1}.Rsquared.Adjusted;
    mlr_roi_model.Yfitted(:,k) = mlr_roi_model.lr_model{k,1}.Fitted;
    mlr_roi_model.residuals_raw(:,k) = mlr_roi_model.lr_model{k,1}.Residuals.Raw;
    mlr_roi_model.residuals_std(:,k) = mlr_roi_model.lr_model{k,1}.Residuals.Standardized;
    mlr_roi_model.residuals_student(:,k) = mlr_roi_model.lr_model{k,1}.Residuals.Studentized;
    mlr_roi_model.ms_error(k,1) = mlr_roi_model.lr_model{k,1}.MSE;
end
mlr_roi_model.X = X;
mlr_roi_model.Y = Y_mi_avg;
%
X1 = X(:,1:3);
Y1 = mlr_roi_model.Yfitted;
clearvars -except X1 Y1 mlr_roi_model
mlr_roi_model_regressed = struct;
for k=1:size(Y1,2)
    fprintf('Fitting the linear regression model for rsn group %d...\n',k);
    mlr_roi_model_regressed.lr_model{k,1} = fitlm(X1,Y1(:,k),'RobustOpts','ols');
    mlr_roi_model_regressed.coeff(k,:) = mlr_roi_model_regressed.lr_model{k,1}.Coefficients.Estimate;
    mlr_roi_model_regressed.pval(k,:) = mlr_roi_model_regressed.lr_model{k,1}.Coefficients.pValue;
    mlr_roi_model_regressed.tstat(k,:) = mlr_roi_model_regressed.lr_model{k,1}.Coefficients.tStat;
    mlr_roi_model_regressed.Rsquared_orig(k,:) = mlr_roi_model_regressed.lr_model{k,1}.Rsquared.Ordinary;
    mlr_roi_model_regressed.Rsquared_adjust(k,:) = mlr_roi_model_regressed.lr_model{k,1}.Rsquared.Adjusted;
    mlr_roi_model_regressed.Yfitted(:,k) = mlr_roi_model_regressed.lr_model{k,1}.Fitted;
    mlr_roi_model_regressed.residuals_raw(:,k) = mlr_roi_model_regressed.lr_model{k,1}.Residuals.Raw;
    mlr_roi_model_regressed.residuals_std(:,k) = mlr_roi_model_regressed.lr_model{k,1}.Residuals.Raw;
    mlr_roi_model_regressed.ms_error(k,1) = mlr_roi_model_regressed.lr_model{k,1}.MSE;
end
mlr_roi_model_regressed.X = X1;
mlr_roi_model_regressed.Y = Y1;
tstat(:,1) = mlr_roi_model.tstat(:,4);
tstat(:,2) = mlr_roi_model_regressed.tstat(:,4);
pval(:,1) = mlr_roi_model.pval(:,4);
pval(:,2) = mlr_roi_model_regressed.pval(:,4);
save(fullfile('E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\Linear',...
    sprintf('mlr_roi_results_sf_zscore.mat')),'mlr_roi_model','mlr_roi_model_regressed','tstat','region_list','pval');
% Plotting the distribution of the residuals
[f1,xi1] = ksdensity(mlr_roi_model_regressed.residuals_raw(:,4));
figure;
plot(xi1,f1)
%lgd = legend('Phonemic fluency');
title('Residuals : Interaction term');
% QQ-plot
figure;
qqplot(mlr_roi_model_regressed.residuals_std(:,4));
title('Residuals : Interaction term');
%------------------------------------------------------------------------------%
%                              Scatter plot - Yfitted vs Fluency score
%------------------------------------------------------------------------------%
% For Left executive control network ROI
cd('E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\Linear');
load('mlr_roi_results_averaged.mat');
figure;
gscatter(mlr_roi_model_regressed.X(:,1),abs(mlr_roi_model_regressed.Y(:,3)),mlr_roi_model_regressed.X(:,2));
title('LECN ROI');
ylabel('Mean FC value');
xlabel('Fluency ratio');
hold on;
plot(mlr_roi_model_regressed.X(1:51,1),mlr_roi_model_regressed.Yfitted(1:51,3),'-r','LineWidth',1);
legend('off');
plot(mlr_roi_model_regressed.X(52:end,1),mlr_roi_model_regressed.Yfitted(52:end,3),'-c','LineWidth',1);
legend('off');
legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
hold off;
% For language ROI
figure;
gscatter(mlr_roi_model_regressed.X(:,1),mlr_roi_model_regressed.Yfitted(:,4),mlr_roi_model_regressed.X(:,2));
title('Language fROI');
ylabel('Mean ROI connectivity');
xlabel('Semantic fluency (z-score)');
hold on;
plot(mlr_roi_model_regressed.X(1:51,1),mlr_roi_model_regressed.Yfitted(1:51,4),'-r','LineWidth',1);
legend('off');
plot(mlr_roi_model_regressed.X(52:end,1),mlr_roi_model_regressed.Yfitted(52:end,4),'-c','LineWidth',1);
legend('off');
legend({'NC', 'MCI'}, 'Location','southoutside');
hold off;
% For precuneus ROI
figure;
gscatter(mlr_roi_model_regressed.X(:,1),mlr_roi_model_regressed.Y(:,8),mlr_roi_model_regressed.X(:,2));
title('Precuneus ROI');
ylabel('Mean FC value');
xlabel('Fluency ratio');
hold on;
plot(mlr_roi_model_regressed.X(1:51,1),mlr_roi_model_regressed.Yfitted(1:51,8),'-r','LineWidth',1);
legend('off');
plot(mlr_roi_model_regressed.X(52:end,1),mlr_roi_model_regressed.Yfitted(52:end,8),'-c','LineWidth',1);
legend('off');
legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
hold off;
% For motor ROI
figure;
gscatter(mlr_roi_model_regressed.X(:,1),mlr_roi_model_regressed.Y(:,5),mlr_roi_model_regressed.X(:,2));
title('Motor ROI');
ylabel('Mean FC value');
xlabel('Fluency ratio');
hold on;
plot(mlr_roi_model_regressed.X(1:51,1),mlr_roi_model_regressed.Yfitted(1:51,5),'-r','LineWidth',1);
legend('off');
plot(mlr_roi_model_regressed.X(52:end,1),mlr_roi_model_regressed.Yfitted(52:end,5),'-c','LineWidth',1);
legend('off');
legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
hold off;
% For primary visual ROI
figure;
gscatter(mlr_roi_model_regressed.X(:,1),mlr_roi_model_regressed.Y(:,10),mlr_roi_model_regressed.X(:,2));
title('dDMN ROI');
ylabel('Mean FC value');
xlabel('Fluency ratio');
hold on;
plot(mlr_roi_model_regressed.X(1:51,1),mlr_roi_model_regressed.Yfitted(1:51,10),'-r','LineWidth',1);
legend('off');
plot(mlr_roi_model_regressed.X(52:end,1),mlr_roi_model_regressed.Yfitted(52:end,10),'-c','LineWidth',1);
legend('off');
legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
hold off;

%------------------------------------------------------------------------------%
%                   Non-linear analysis : Randomforest regression
%------------------------------------------------------------------------------%
% Generating the cross validation partition

Data = array2table(cv.X_train);
Data.Yfit = cv.Y_train(:,1);
t = templateTree('NumVariablesToSample','all',...
    'PredictorSelection','interaction-curvature','Surrogate','on');
rng(1); % For reproducibility
Mdl = fitrensemble(X,Y_mean(:,1),'Method','Bag','NumLearningCycles',200, ...
    'Learners',t);

load imports-85
Y = X(:,1);
X = X(:,2:end);
cat_var = [0 1 0 0 0 0 0];
b = TreeBagger(200,X,Y_mean(:,1),'Method','regression', ...
    'OOBPredictorImportance','On', ...
    'CategoricalPredictors',find(cat_var == 1), ...
    'MinLeafSize',100);

Mdl = TreeBagger(30, Data, 'Yfit', ...
    'Method', 'regression', 'PredictorSelection', 'curvature', ...
    'Surrogate', 'on', 'OOBPredictorImportance', 'on');

%Trial - 
% The entire time series
r_pearson = corrcoef(timeseries_region);
r_spearmann = corr(timeseries_region,'type','Spearman', 'rows','complete');
r_kendall = corr(timeseries_region, 'type','Kendall');






