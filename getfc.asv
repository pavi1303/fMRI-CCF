function FCN_VEC = getfc(fcn_mat,idx_start,idx_end)
% Making the upper triangular matrix NaN
temp=ones(size(fcn_mat));
%idx=tril(temp,-1);
%fcn_mat(~idx)=nan;
fcn_mat(eye(size(fcn_mat))==1) = nan;
% Extracting the corresponding columns without the NaN values
fcn_col = fcn_mat(:,idx_start:idx_end);
FCN_COL = reshape(fcn_col.',1,[]);
% Returning the vector of functional connections for that rsn group
FCN_VEC = FCN_COL(~isnan(FCN_COL));
end