function [fcn_mean, ] = fcnmat_results(fcn_dir, pos_thresh, neg_thresh, n_ica)
dirloc = dir(fcn_dir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i = 1:length(subloc)
    sub = subloc{i};
    current = strcat(fcn_savedir,'\',sub);
    sub_data = load(current,'fcn');
    fcn_mat{1,i} = (sub_data.fcn);
end
X = cat(3,fcn_mat{:});
fcn_mean = mean(X,3);
L_mat = tril(fcn_mean,-1);
[row, col, ~] = find(L_mat>pos_thresh | L_mat<neg_thresh);
idx = horzcat(row, col);
val = zeros(n_ica,n_ica);
for i = 1:size(idx,1)
    val(idx(i,1),idx(i,2)) = L_mat(idx(i,1),idx(i,2));
end

end
