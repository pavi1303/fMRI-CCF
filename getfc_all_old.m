function [fcn_vec, idx] = getfc_all(fcn_mat)
% Making the upper triangular matrix NaN
temp=ones(size(fcn_mat));
idx=tril(temp,-1);
fcn_mat(~idx)=nan;
fcn_mat_nonan = ~isnan(fcn_mat);
[row,col] = find(fcn_mat_nonan);
idx = horzcat(row,col);
fcn_vec = fcn_mat(fcn_mat_nonan(:));
end