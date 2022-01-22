function [data] = generate_vt(list_of_files, mask, xres, yres, zres)
for j=1:length(list_of_files)
    i = load_untouch_nii(list_of_files(j).name);
    I = i.img;
    I_M = I.*mask;
    vt = reshape(I_M, [1,xres*yres*zres]);
    [~,idx] = find(vt);
    vt = vt(:,idx);
    data(j,:) = vt;
end
end