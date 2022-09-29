function fcn_vec = getfcn_vec(A)
mask = tril(true(size(A)),-1);
fcn_vec = A(mask)';
end
