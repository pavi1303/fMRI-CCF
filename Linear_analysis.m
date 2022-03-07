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

% GENERATION OF THE RESPONSE VARIABLE
fcn_savedir='E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\4.FCN_reordered';
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
    fcn_vis = getfc(sub_fcn,1,10);
    fcn_aud = getfc(sub_fcn,11,14);
    fcn_lang = getfc(sub_fcn,15,28);
    fcn_mem = getfc(sub_fcn,29,44);
    fcn_sub = getfc(sub_fcn,45,53);
    fcn_cer = getfc(sub_fcn,54,64);
    fcn_mot = getfc(sub_fcn,65,70);
    fcn_sen = getfc(sub_fcn,71,74);
    fcn_all = horzcat(fcn_vis,fcn_aud,fcn_lang,fcn_mem,fcn_sub,fcn_cer,fcn_mot,fcn_sen);
    fcn_all_mean = horzcat(mean(fcn_vis),mean(fcn_aud),mean(fcn_lang),mean(fcn_mem),...
        mean(fcn_sub),mean(fcn_cer),mean(fcn_mot),mean(fcn_sen));
    Y(j,:) = fcn_all;
    Y_mean(j,:) = fcn_all_mean;
end

clearvars -except X Y Y_mean;
%Performing linear regression by myself
X1 = horzcat(ones(size(X,1),1),X);
XT = X1';
XM = XT*X1;
XINV = inv(XM);
b = (XINV*(XT*Y_mean))';
% Peforming multiple linear regression model - averaged across groups
for k=1:size(Y_mean,2)
    fprintf('Fitting the linear regression model for rsn group %d...\n',k);
    lr_model{k,1} = fitlm(X,Y_mean(:,k),'RobustOpts','ols');
    coeff(k,:) = lr_model{k,1}.Coefficients.Estimate;
    pval(k,:) = lr_model{k,1}.Coefficients.pValue;
    tstat(k,:) = lr_model{k,1}.Coefficients.tStat;
    Rsquared_orig(k,:) = lr_model{k,1}.Rsquared.Ordinary;
    Rsquared_adjust(k,:) = lr_model{k,1}.Rsquared.Adjusted;
    Yfitted(:,k) = lr_model{k,1}.Fitted;
    residuals_raw(:,k) = lr_model{k,1}.Residuals.Raw;
    residuals_std(:,k) = lr_model{k,1}.Residuals.Raw;
    ms_error(k,1) = lr_model{k,1}.MSE;
end
% Fitting multiple linear regression model - individual connections
for k=1:size(Y,2)
    fprintf('Fitting the linear regression model for rsn group %d...\n',k);
    lr_model{k,1} = fitlm(X,Y(:,k),'RobustOpts','ols');
    coeff(k,:) = lr_model{k,1}.Coefficients.Estimate;
    pval(k,:) = lr_model{k,1}.Coefficients.pValue;
    tstat(k,:) = lr_model{k,1}.Coefficients.tStat;
    Rsquared_orig(k,:) = lr_model{k,1}.Rsquared.Ordinary;
    Rsquared_adjust(k,:) = lr_model{k,1}.Rsquared.Adjusted;
    Yfitted(:,k) = lr_model{k,1}.Fitted;
    residuals_raw(:,k) = lr_model{k,1}.Residuals.Raw;
    residuals_std(:,k) = lr_model{k,1}.Residuals.Raw;
    ms_error(k,1) = lr_model{k,1}.MSE;
end
avg_tstat(1,1) = mean(tstat(1:685,4));
avg_tstat(1,2) = mean(tstat(686:931,4));
avg_tstat(1,3) = mean(tstat(932:1666,4));
avg_tstat(1,4) = mean(tstat(1667:2266,4));
avg_tstat(1,5) = mean(tstat(2267:2491,4));
avg_tstat(1,6) = mean(tstat(2492:2656,4));
avg_tstat(1,7) = mean(tstat(2657:2695,4));
avg_tstat(1,8) = mean(tstat(2696:2701,4));
%------------------------------------------------------------------------------%
%             BOOTSTRAP AGGREGATION - USING TREBAGGER
%------------------------------------------------------------------------------%
% 1. Finding the optimal leaf size
for j = 1:size(Y_mean,2)
    Y1 = Y_mean(:,j);
    isCategorical = [0,1,0,0,0,0];
    leaf = [5 10 20 50 100];
    col = 'rbcmk';
    figure
    hold on
    for i=1:length(leaf)
        b = TreeBagger(50,X,Y1,'Method','regression', ...
            'OOBPrediction','On', ...
            'CategoricalPredictors',find(isCategorical == 1), ...
            'MinLeafSize',leaf(i));
        plot(oobError(b),col(i))
    end
    title('Network %d',j);
    xlabel('Number of Grown Trees');
    ylabel('Mean Squared Error');
    legend({'5' '10' '20' '50' '100'},'Location','NorthEast');
    hold off
end

% 2. Estimating feature importance
for j = 1:size(Y_mean,2)
    b = TreeBagger(200,X,Y_mean(:,j),'Method','regression', ...
        'OOBPredictorImportance','On', ...
        'CategoricalPredictors',find(isCategorical == 1), ...
        'MinLeafSize',100);
    figure;
    plot(oobError(b));
    xlabel('Number of Grown Trees');
    ylabel('Out-of-Bag Mean Squared Error');
end


b1 = fitrensemble(X,Y(:,1),'Method','Bag');
figure;
bar(b.OOBPermutedPredictorDeltaError);
xlabel('Feature Number') ;
ylabel('Out-of-Bag Feature Importance');

Imp = oobPermutedPredictorImportance(b);

for k=1:size(Y,2)
    fprintf('Fitting the non-linear regression model for rsn group %d...\n',k);
    rf_model{k,1} = TreeBagger(100,X,Y(:,k),'Method','regression');
end

X = rand(106,6);
Y_mean = rand(106,8);
rng(1945,'twister')
for k=1:size(Y_mean,2)
    fprintf('Fitting the randomforest regression model for rsn group %d...\n',k);
    rf_model{k,1} = TreeBagger(50,X,Y_mean(:,k),'Method','regression','PredictorSelection','curvature',...
        'OOBPredictorImportance','on');
end

%%%%%%%%%%% COMPARING THE FLUENCY RATIO
load carsmall
Cylinders = categorical(Cylinders);
Mfg = categorical(cellstr(Mfg));
Model_Year = categorical(Model_Year);
X = table(Acceleration,Cylinders,Displacement,Horsepower,Mfg,...
    Model_Year,Weight,MPG);
rng('default'); % For reproducibility
Mdl = TreeBagger(200,X,'MPG','Method','regression','Surrogate','on',...
    'PredictorSelection','curvature','OOBPredictorImportance','on');