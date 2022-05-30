function [fcn_values, fcn_vec] = getfc_all(fcn_mat)
% Extracting all the values
fcn_ltri_mat= tril(fcn_mat, -1);
fcn_values = (abs(fcn_ltri_mat(fcn_ltri_mat~=0)))';
% Making the upper triangular matrix NaN
temp=ones(size(fcn_mat));
fcn_mat(eye(size(fcn_mat))==1) = nan;
% Extracting the corresponding columns without the NaN values
for col_idx = 1:size(fcn_mat,2)
    fcn_col = fcn_mat(:,col_idx);
    fcn_vec(1,col_idx) = mean(fcn_col(~isnan(fcn_col)));
    fcn_vec = abs(fcn_vec);
end
end
