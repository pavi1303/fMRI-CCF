function [data_red] = do_PCA(data, n_comp)
    c = cov(data');
    [eigvec,eigenval,~]=pcacov(c);
    EV=eigvec(:,1:n_comp);
    data_red=EV'*data;
end