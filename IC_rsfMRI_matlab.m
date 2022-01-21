%Getting the required subdir
clc
clear
tic
rootdir = 'E:\LRCBH\Data\COBRE-MNI\Trial';
savedir = 'E:\LRCBH\Results\Matlab\1.PCA';
dirpath = dir(rootdir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
%Loading the mask file
cd('E:\LRCBH\MNI_segmented');
m = load_untouch_nii('standard_binary.nii');
M = m.img;
[x, y, z] = size(M);
trs = 850;
comp = 200;
%Iterating through each of the subjects
for i=1:length(subloc)
    suboi = subloc{i};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    files = dir('*.nii');
    temp = zeros(trs, 228453);
    fprintf('Generating voxeltime data for subject %s...\n',suboi);
    temp = double(generate_vt(files,M,x,y,z));
    vt_data = zscore(temp(16:end, :),1);
    clear temp;
    fprintf('Performing PCA reduction using %d components for subject %s...\n',comp,suboi);
    %Performing subject wise PCA reduction
    PCA_red = double(do_PCA(vt_data,comp));
    cd(savedir);
    save(fullfile(savedir, sprintf('PCA_%s.mat',suboi)),'PCA_red');
end
fprintf('PCA reduction done.\n');
toc
