%------------------------------------------------------------------------------%
%         DATA PREPARATION - WILLARD ATLAS BASED ANALYSIS
%------------------------------------------------------------------------------%

rootdir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data\Useful';
atlas_loc = 'E:\LRCBH\Atlas\Willard_with_overlap_2mm\1.Merged';
data_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\1.Data_merged';
fcn_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\3.FCN\1.Merged';

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
            [~,idx] = find(vt);
            vt = vt(:,idx);
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
interaction = readmatrix('Cobre_fluency_study_v2.xlsx', 'Sheet','regression','Range',[2 5 107 5]);
%suvr = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 8 103 8]);
%pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 9 103 9]);
%sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 10 103 10]);
suvr_dis = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 11 107 11]);

%regressor = horzcat(fluency_ratio, grp);
regressor = horzcat(fluency_ratio, grp, interaction);
covariates = horzcat(age, ed, suvr_dis);
X = horzcat(regressor, covariates);

% GENERATION OF THE FUNCTIONAL CONNECTIVITY MATRIX
dirpath = dir(data_savedir);
subpath = {dirpath.name}';
subpath(ismember(subpath,{'.','..'})) = [];
for i=1:length(subpath)
    sub = subpath{i};
    current = strcat(data_savedir,'\', sub);
    tc_data = load(current,'timeseries_region');
    tc = (tc_data.timeseries_region);
    [pval, fcn_roi] = corrcoef(tc);
    save(fullfile(fcn_savedir, sprintf('FCN_merged_%s.mat',sub)),'fcn_roi','pval');
end
% GENERATION OF THE RESPONSE VARIABLE
dirloc = dir(fcn_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for j = 1:length(subloc)
    sub = subloc{j};
    current = strcat(fcn_savedir,'\',sub);
    sub_fcn = load(current,'fcn_roi');
    sub_fcn = sub_fcn.fcn_roi;
    fcn_aud = getfc_v1(sub_fcn,1);
    fcn_bas = getfc_v1(sub_fcn,2);
    fcn_lecn = getfc_v1(sub_fcn,3);
    fcn_lang = getfc_v1(sub_fcn,4);
    fcn_mot = getfc_v1(sub_fcn,5);
    fcn_prec = getfc_v1(sub_fcn,6);
    fcn_recn = getfc_v1(sub_fcn,7);
    fcn_sal = getfc_v1(sub_fcn,8);
    fcn_visspa = getfc_v1(sub_fcn,9);
    fcn_ddmn = getfc_v1(sub_fcn,10);
    fcn_hvis = getfc_v1(sub_fcn,11);
    fcn_psal = getfc_v1(sub_fcn,12);
    fcn_pvis = getfc_v1(sub_fcn,13);
    fcn_vdmn = getfc_v1(sub_fcn,14);
    fcn_all = horzcat(fcn_aud,fcn_aud,fcn_bas,fcn_lecn,fcn_lang,fcn_mot,fcn_mot,fcn_prec,...
       fcn_recn,fcn_sal,fcn_visspa,fcn_ddmn,fcn_hvis,fcn_psal,fcn_pvis,fcn_vdmn);
    fcn_all_mean = horzcat(mean(fcn_aud),mean(fcn_bas),mean(fcn_lecn),...
        mean(fcn_lang),mean(fcn_mot),mean(fcn_prec),mean(fcn_recn),...
        mean(fcn_sal),mean(fcn_visspa),mean(fcn_ddmn),mean(fcn_hvis),mean(fcn_psal),...
        mean(fcn_pvis),mean(fcn_vdmn));
    Y(j,:) = fcn_all;
    Y_mean(j,:) = fcn_all_mean;
end
clearvars -except X Y Y_mean ;
% Fitting the regression model
mlr_roi_model= struct;
for k=1:size(Y_mean,2)
    fprintf('Fitting the linear regression model for rsn group %d...\n',k);
    mlr_roi_model.lr_model{k,1} = fitlm(X,Y_mean(:,k),'RobustOpts','ols');
    mlr_roi_model.coeff(k,:) = mlr_roi_model.lr_model{k,1}.Coefficients.Estimate;
    mlr_roi_model.pval(k,:) = mlr_roi_model.lr_model{k,1}.Coefficients.pValue;
    mlr_roi_model.tstat(k,:) = mlr_roi_model.lr_model{k,1}.Coefficients.tStat;
    mlr_roi_model.Rsquared_orig(k,:) = mlr_roi_model.lr_model{k,1}.Rsquared.Ordinary;
    mlr_roi_model.Rsquared_adjust(k,:) = mlr_roi_model.lr_model{k,1}.Rsquared.Adjusted;
    mlr_roi_model.Yfitted(:,k) = mlr_roi_model.lr_model{k,1}.Fitted;
    mlr_roi_model.residuals_raw(:,k) = mlr_roi_model.lr_model{k,1}.Residuals.Raw;
    mlr_roi_model.residuals_std(:,k) = mlr_roi_model.lr_model{k,1}.Residuals.Raw;
    mlr_roi_model.ms_error(k,1) = mlr_roi_model.lr_model{k,1}.MSE;
end
mlr_roi_model.X = X;
mlr_roi_model.Y = Y_mean;
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
    sprintf('mlr_roi_results_averaged.mat')),'mlr_roi_model','mlr_roi_model_regressed','tstat','region_list','pval');
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


