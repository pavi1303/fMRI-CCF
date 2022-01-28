function [X,Y,lr_model,coeff,pval,tstat,alpha_level,pval_interaction,Rsquared_orig_mean,Rsquared_adjust_mean,sig_voxel] = regress_model(dr_dir, regressor, interaction, covariates, var_name,rsn_no,save_prefix,savepath)
% Initializing the alpha value
alpha_level = 0.05;
% Loading the subject specific spatial maps
dirloc = dir(dr_dir);
subdir = [dirloc(:).isdir];
subloc = {dirloc(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(dr_dir,'\', suboi);
    cd(current);
    sub_sm = load(var_name);
    sub_sm = (sub_sm.ss_sm)';
    sub_smoi(i,:) = sub_sm(rsn_no,:);
end
% Forming the Y data
Y = sub_smoi;
%Forming the data matix
X = horzcat(regressor,interaction,covariates);
% Preforming regression
for k=1:size(Y,2)
    lr_model{k,1} = fitlm(X,Y(:,k),'RobustOpts','ols');
    coeff(k,:) = lr_model{k,1}.Coefficients.Estimate;
    pval(k,:) = lr_model{k,1}.Coefficients.pValue;
    tstat(k,:) = lr_model{k,1}.Coefficients.tStat;
    Rsquared_orig(k,:) = lr_model{k,1}.Rsquared.Ordinary;
    Rsquared_adjust(k,:) = lr_model{k,1}.Rsquared.Adjusted;
end
%Finding the mean coefficient of determination
Rsquared_orig_mean = mean(Rsquared_orig);
Rsquared_adjust_mean = mean(Rsquared_adjust);
% Correcting for multiple comparisons
% n_compar = size(Y,2);
pval_interaction = (pval(:,4))';% For the interaction term
%pval_measure = (pval(:,1))';% For the fluency term
%alpha_level = 0.05/n_compar;
sig_idx= find(pval_interaction<alpha_level);
sig_voxel = vertcat(sig_idx,pval_interaction(:,sig_idx));
save(fullfile(savepath, sprintf('regress_model_%s.mat',save_prefix)),'X','Y','lr_model','coeff','pval','tstat','Rsquared_orig','Rsquared_adjust','Rsquared_orig_mean','Rsquared_adjust_mean','alpha_level','pval_interaction','sig_voxel');
end
