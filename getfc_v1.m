function FCN_VEC = getfc_v1(fcn_mat,col_idx)
% Making the upper triangular matrix NaN
temp=ones(size(fcn_mat));
fcn_mat(eye(size(fcn_mat))==1) = nan;
% Extracting the corresponding columns without the NaN values
fcn_col = fcn_mat(:,col_idx);
FCN_COL = reshape(fcn_col.',1,[]);
% Returning the vector of functional connections for that rsn group
FCN_VEC = abs(FCN_COL(~isnan(FCN_COL)));
end