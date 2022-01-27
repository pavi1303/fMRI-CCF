function [Yfitted, Yfitted_sig,beta,pval,pval_sig,tstat,alpha_level,sig_loc] = confound_fitlm(fcn_dir, X, var_name,no_rsn,save_prefix,savepath)
% Initializing the alpha value
alpha_level = 0.05;
% Finding the non-zero networks
c = corrcoef(rand(no_rsn,835)');
corr_mat = tril(c,-1);
[row,col] = find(corr_mat);
nz_loc = horzcat(row,col);
% Loading the functional connectivity matrix
dirloc = dir(fcn_dir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(fcn_dir,'\', suboi);
    sub_data = load(current,var_name);
    sub_data = (sub_data.fcn)';
    sub_data(2,:) = [];
    sub_data(:,2) = [];
    corr_mat = tril(sub_data,-1);
    corr_vec = nonzeros(corr_mat);
    fnc_val{i,1} = corr_vec';
end
% Forming the Y data
Y = vertcat(fnc_val{:});
% Preforming regression
for k=1:size(Y,2)
    lr_model{k,1} = fitlm(X,Y(:,k),'RobustOpts','ols');
    beta(k,:) = lr_model{k,1}.Coefficients.Estimate;
    pval(k,:) = lr_model{k,1}.Coefficients.pValue;
    Yfitted(:,k) = lr_model{k,1}.Fitted;                                 
    tstat(k,:) = lr_model{k,1}.Coefficients.tStat;
end
% Correcting for multiple comparisons
n_compar = size(Y,2);
p_val_measure = (pval(:,2))';
%alpha_level = 0.05/n_compar;
sig_idx = find(p_val_measure<alpha_level);
pval_sig = p_val_measure(:,sig_idx);
Yfitted_sig = Yfitted(:,sig_idx);
sig_loc = nz_loc(sig_idx,:);
save(fullfile(savepath, sprintf('fcn_association_%s.mat',save_prefix)),'beta','pval','tstat','alpha_level','sig_loc');
end
