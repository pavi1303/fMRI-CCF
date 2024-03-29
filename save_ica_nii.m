function save_ica_nii(ica_mat,xdim,ydim,zdim,idx,nii_temp,save_prefix,savepath)
if ~exist(savepath,'dir')
    mkdir(savepath);
end
for j = 1:size(ica_mat,1)
    temp = zeros(1, xdim*ydim*zdim);
    ica_comp = ica_mat(j,:);
    temp(1,idx) = ica_comp;
    ica_3d = reshape(temp,[xdim,ydim,zdim]);
    nii_temp.img = [];
    nii_temp.img = single(ica_3d);
    
    save_untouch_nii(nii_temp,[savepath,'\',Sname]);
end