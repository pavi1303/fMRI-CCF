function pcorr = partialcor(data)
pres_mat = inv(cov(data));
for i=1:size(pres_mat)
    for j=1:size(pres_mat)
        pcorr(i,j) = -(pres_mat(i,j)/sqrt(pres_mat(i,i)*pres_mat(j,j)));
    end
end
pcorr(eye(size(pcorr))==1) = 1;
end

