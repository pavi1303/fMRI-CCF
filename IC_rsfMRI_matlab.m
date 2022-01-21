root_dir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data';
subdir = [root_dir(:).isdir];
nameFolds = {root_dir(subdir).name}';
subfolder = string(nameFolds);
nameFolds(ismember(nameFolds,{'.','..'})) = [];
files = dir(fullfile(strcat(root_dir,'\',subfolder(1,1))), '*.nii');

addpath(genpath('E:\LRCBH\Data\COBRE-MNI\Trial\008'))
addpath(genpath('E:\LRCBH\MNI_segmented'))

data = load_untouch_nii('s8_MNI_012_00001.nii');
img_data = data.img;
rootdir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data';
filelist = dir(fullfile(rootdir, '**\*.*'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]); 

cd('C:\Users\PATTIAP\Desktop\Dataset')
data = load_untouch_nii('MNI-008_unmasked.nii');
mask = load_untouch_nii('standard_binary.nii');
mask_data = mask.img;
vt_data = xyz_to_q_ALT(data.img);
%Applying mask
masked_data = bsxfun(@times, img_data, mask_data);
[x,y,z] = size(masked_data);
masked_img = permute(masked_data, [4, 1, 2, 3]);
vt_data = reshape(masked_data, [1, x*y*z]);
vt_data_modified = nonzeros(vt_data)';
ig_size = masked_img.
%%% For iterating through different folders
basePath = pwd;  %your base path which is in your case myTraining  
allPaths = dir(basePath);  %get all directory content
subFolders = [allPaths(:).isdir]; %get only indices of folders
foldersNames = {allPaths(subFolders).name}'; % filter folders names
foldersNames(ismember(foldersNames,{'.','..'})) = []; %delete default paths for parents return '.','..'
for i=1:length(foldersNames), %loop through all folders
    tmp = foldersNames{i};  %get folder by index
    p = strcat([basePath '\']); 
    currentPath =strcat([p tmp]); % add base to current folder
    cd(currentPath);   % change directory to new path
    files = dir('*.jpg'); % list all images in your path which in your case could be John or Mary 
    for j=1:length(files), % loop through your images 
        img = imread(files(j).name); % read each image and do what you want 
    end
end 