function [fluency_ratio, subloc, fcn_mat, rsn_idx, corr_value, corr_min, corr_max,N, poscorr_value,negcorr_value,rsn_pos,rsn_neg] = ...
    fcnmat_results(fcn_dir, fluency_loc, pat_idx, noise_idx,grp_id, xrange, poscorr_range, negcorr_range)
% Loading the fluency ratio 
cd(fluency_loc);
if grp_id ==1
    fluency_ratio = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 12 52 12]);
elseif grp_id==2
    fluency_ratio = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[53 12 103 12]);
elseif grp_id==3
    fluency_ratio = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 12 103 12]);
else
    fprintf('Invalid choice\n');
end
% Getting the list of patients in the corresponding group
dirloc = dir(fcn_dir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
subloc = subloc(pat_idx,:);
% Generating the patient * fcn value matrix 
for i = 1:length(subloc)
    sub = subloc{i};
    current = strcat(fcn_dir,'\',sub);
    sub_fcn = load(current,'fcn');
    sub_fcn = sub_fcn.fcn;
    sub_fcn(:,noise_idx) = 0;
    sub_fcn(noise_idx,:) = 0;
    L_mat = tril(sub_fcn,-1);
    fcn_mat(i,:) = L_mat(L_mat~=0);
end
% Getting the indices of the resting state networks
[row,col] = find(L_mat);
rsn_idx = horzcat(row,col);
% Getting the correlation values
for j = 1:size(fcn_mat,2)
    fcn_oi = fcn_mat(:,j);
    r = corrcoef(fcn_oi,fluency_ratio);
    corr_value(1,j) = r(1,2);
end
% Finding and returning the max and min corr values
corr_value = (corr_value)';
corr_min = min(corr_value);
corr_max = max(corr_value);
% Finding the connections in a specific range and it's value
N = histcounts(corr_value, xrange);
pos_idx = find(poscorr_range(1,1)<=corr_value & corr_value<poscorr_range(1,2));
neg_idx = find(negcorr_range(1,1)<=corr_value & corr_value<negcorr_range(1,2));
poscorr_value = corr_value(pos_idx,1);
negcorr_value = corr_value(neg_idx,1);
% Finding the networks with those values
rsn_pos = rsn_idx(pos_idx,:);
rsn_neg = rsn_idx(neg_idx,:);
end