%Getting the required subdir
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
%Iterating through each of the subjects
for i=1:length(subloc)
    suboi = subloc{i};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    files = dir('*.nii');
    temp = zeros(trs, 228453);
    for j=1:length(files)
        i = load_untouch_nii(files(j).name);
        I = i.img;
        I_M = I.*M;
        vt = reshape(I_M, [1,x*y*z]);
        [~,idx] = find(vt);
        vt = vt(:,idx);
        temp(j,:) = vt;
        clear vt, i, I;
    end
    vt_data = zscore(temp(16:end, :));
    clear temp;
    %Performing subject wise PCA reduction
    vt_data= vt_data - mean(vt_data,2);
    c = cov(vt_data');
    [V,D,explained]=pcacov(c);
    EV=V(:,1:200);
    PCA_red=EV'*vt_data;
    cd(savedir);
    save(fullfile(savedir, sprintf('PCA_%s.mat',suboi)),'PCA_red');
end

i=1;

files = dir('*.nii');
temp = zeros(trs, 228453);
for j = 1:length(files)
    i = load_untouch_nii(files(j).name);
    I = i.img;
    I_M = I.*M;
    vt = reshape(I_M, [1,x*y*z]);
    [~,idx] = find(vt);
    vt = vt(:,idx);
    temp(j,:) = vt;
    clear vt, i, I;
end
vt_data = zscore(temp(16:end, :));
clear temp
%Performing subject wise PCA reduction
vt_data= vt_data - mean(vt_data,2);
c = cov(vt_data');
[V,D,explained]=pcacov(c);
EV=V(:,1:200);
PCA_red=EV'*vt_data;
cd(savedir);
savename = strcat()
save(strcat(savedir,'\','PCA%s.mat',suboi),'PCA_red)';

savename = fullfile(savedir, sprintf('PCA%s.mat',suboi));
savename = fullfile(savedir, sprintf('PCA_%s.mat',suboi));
save(fullfile(savedir, sprintf('PCA_%s.mat',suboi)),'PCA_red');

baseFileName = sprintf('%s.mat', sixCenters{k});
fullMatFileName = fullfile(folder, baseFileName);
save(fullMatFileName, 'grouped');


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
masked_data = bsxfun(@times, img, mask_data);
[x,y,z] = size(mask_data);
masked_img = permute(masked_data, [4, 1, 2, 3]);
vt_data = reshape(masked_data, [1, x*y*z]);
vt_data_modified = nonzeros(vt_data)';
vt_data_modified = zscore(vt_data_modified);
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

% My try 1 for masking the image from a patient subfolder
% Loading of the mask file way outside the for loop and getting its
% dimensions
cd('C:\Users\pavig\Downloads');
m = load_untouch_nii('standard_binary.nii');
M = m.img;
[xres, yres, zres] = size(M);
% Get the number of time series from the user
trs = 850;
% This is the inner for loop where the generation of the voxel time data
% happens and the conversion of the PCA also happens for each subject and
% then I save th result of PCA as a mat variable. 
cd('W:\008');
files = dir('*.nii');
temp = zeros(trs, xres*yres*zres);
for j = 1:length(files)
    i = load_untouch_nii(files(j).name);
    I = i.img;
    I_M = I.*M;
    vt = reshape(I_M, [1,xres*yres*zres]);
    %vt = nonzeros(vt)';
    [~,idx] = find(vt); 
    temp(j,:) = vt;
    clear vt, i, I;
end
idx = unique(v);
[~,v] = find(temp);
vt_data = temp(16:end, :);
X = vt_data;
%PCA analysis in matlab
X = X - mean(X,2);
c = cov(X');
[E,eigval,explained]=pcacov(c); 
EE=E(:,1:200);
X_tilde=EE'*X;

mean_x=mean(X,2);
[coeff,score] = pca(vt_data);
reduced = score(
    
    
    
% You load the mask image outside the for loop
j =1;
data = load_untouch_nii(files(j).name);
img = data.img;
img_masked = img.*mask_data;
vt_data = reshape(I_M, [1, xres*yres*zres]);
vt = nonzeros(vt_data)';
vt = zscore(vt);

% Items to do
%[1] Load the ydata file from xz and make sure that the PCA process is being done
%correctly and also check the mask computation process and the indexes for
%it. 
%[2] Include to save the voxeltime data in the file and then see if the
%indexes change - they shouldn't coz they are from the mask whose 0 and 1 locations are the same.
%[3] Use these indexes to then reconstruct the 4D array from the 2D data. 
%[4] After you figure out the reconstruction part, then  prepare the code
%getting the individual slices of the 4D ICA data.
% [5] Test this out on the trial dataset. If it works - then now run the
% program and go to sleep. Check the results tomorrow. 