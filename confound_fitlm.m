function regress_model(dr_dir, regressor, interaction, covariates, var_name,rsn_no,save_prefix,savepath)
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
X = horzcat(regressor,interaction,covariates);% With interaction
% Preforming regression
for k=1:size(Y,2)
    lr_model{k,1} = fitlm(X,Y(:,k),'RobustOpts','ols');
    beta(k,:) = lr_model{k,1}.Coefficients.Estimate;
    pval(k,:) = lr_model{k,1}.Coefficients.pValue;            
    tstat(k,:) = lr_model{k,1}.Coefficients.tStat;
end
% Correcting for multiple comparisons
% n_compar = size(Y,2);
pval_interaction = (pval(:,4))';% For the interaction term
%pval_measure = (pval(:,1))';% For the fluency term
%alpha_level = 0.05/n_compar;
sig_idx= find(p_val_interaction<alpha_level);
sig_voxel = vertcat(sig_idx,p_val_interaction(:,sig_idx));
save(fullfile(savepath, sprintf('sm_association_%s.mat',save_prefix)),'lr_model','beta','pval','tstat','alpha_level','pval_interaction','sig_voxel');
end
