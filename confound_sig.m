function [coeff,pval_corr,stats,n_compar,sig_loc] = confound_sig(fcn_dir, X, var_name,corrmat_template)
% Finding the non-zero networks
corr_mat = tril(corrmat_template,-1);
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
    corr_mat = tril(sub_data,-1);
    corr_vec = nonzeros(corr_mat);
    fnc_val{i,1} = corr_vec';
end
% Forming the Y data
Y = vertcat(fnc_val{:});
% Preforming regression
for k = 1:size(Y,2)
    [coeff(:,k),~,~,~,stats(:,k)] = regress(Y(:,k),X);
end
% Correcting for multiple comparisons
n_compar = size(Y,2);
stats(3,:) = stats(3,:)/n_compar;
pval_corr = 0.05/n_compar;
sig_idx = find(stats(3,:)<pval_corr);
sig_loc = nz_loc(sig_idx,:);
end
