
current = 'E:\LRCBH\Data\COBRE-MNI\Individual_data\Useful\008'
for j = 1:length(subloc)
    suboi = subloc{j};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    filelist = dir('*.nii');
    timeseries_voxel = [];
    timeseries_region = [];
    fprintf('Generating data for subject %s...\n',suboi);
    for i = 1:length(region_list)
        cd(atlas_loc);
        region = load_untouch_nii(region_list{i,1});
        roi = single(region.img);
        [x, y, z] = size(roi);
        roi_temp = reshape(roi, [1,x*y*z]);
        [~, indices] = find(roi_temp);
        for k=1:length(filelist)
            cd(current);
            vol = load_untouch_nii(filelist(k).name);
            I = vol.img;
            I_R = I.*roi;
            vt = reshape(I_R, [1,x*y*z]);
            %[~,idx] = find(vt);
            vt = vt(:,indices);
            data(k,:) = vt;
            avg_data = mean(data,2);
        end
        timeseries_voxel = double(horzcat(timeseries_voxel,data));
        timeseries_region = double(horzcat(timeseries_region,avg_data));
        voxel_loc{i,1} = indices;
        clear vol I I_R vt data avg_data region roi x y z roi_temp indices;
    end
    save(fullfile(data_savedir, sprintf('ROI_%s.mat',suboi)),'timeseries_voxel','timeseries_region','voxel_loc');
    clear timeseries_region timeseries_voxel voxel_loc;
end

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
    Sname = [save_prefix,sprintf(['%0d.nii'],j)];
    save_untouch_nii(nii_temp,[savepath,'\',Sname]);
end
save_prefix = 'gICA'
j=1;
Sname = [save_prefix,sprintf(['%0d.nii'],j)];
% Verifying the time series extraction protocol
mask_loc = 'E:\LRCBH\Data\MNI_segmented';
cd(mask_loc);
mask = load_untouch_nii('GMWM_mask.nii')
M = single(mask.img);
pat_nii = 'E:\LRCBH\Data\COBRE-MNI\Individual_data\Useful\008\';
cd(pat_nii);
filelist = dir('*.nii');
filelist = filelist(16:end, 1);
fmri_img = load_untouch_nii(filelist(1).name);
I = fmri_img.img;
I_M = I.*M;
fmri_temp = fmri_img;
fmri_temp.img = I_M;
save_untouch_nii(fmri_temp,[pat_nii,'\',Sname]);
[x, y, z] = size(M);
vt_M = reshape(M, [1,x*y*z]);
c = unique(vt_M);
[x, y, z] = size(I_M);
vt_I_M = reshape(I_M, [1,x*y*z]);
c = unique(vt_I_M);

L_ROI = 'E:\LRCBH\Atlas\Willard_with_overlap_2mm\1.Merged'
cd(L_ROI);
filelist = dir('*.nii');
lan_roi = load_untouch_nii(filelist(1).name);
% 


ts_normalised = zscore(timeseries_region');

ts1 = timeseries_region(:, 1)';
ts1_n = zscore(ts1);
tsn1 = ts_normalised(1, :);
mu = mean(ts1_n);
sigma = std(ts1_n);
for i = 1: length(ts1_n)
    if tsn1(1, i) < u
        zu(1, i) = 0;
    else
        zu(1, i) = 1;
    end
end

pd = makedist('Normal');
u = cdf(pd, 0.75)
pd = makedist('Normal');
l = cdf(pd, -0.75)
trs = timeseries_region';
% -------------------------------------------------------------------------------- %
%           ACCORDANCE AND DISCORDANCE BETWEEN 2 TIME SERIES
% -------------------------------------------------------------------------------- %
% Quantile q is chosen as 0.75 
% Estimation of the lower and upper limits u and l
% STEP 1: Normalization of the time series
for j = 1:size(trs, 1)
    trs_n(j, :) = zscore(trs(j, :));
end
% STEP 2: dentifying the extreme sub intervals in the entire time series 
% Estimating the upper and lower threshold levels ( u and l) by using the
% CDF of the gaussian distribution
q = 0.75
u = norminv(q);
l = norminv(1-q);
% STEP 3: Forming the activation vector for every time series 
for i = 1:size(trs_n, 1)
    trs = trs_n(i, :)
    for j = 1: length(trs)
        if trs(1, j)<u
            xu(i, j) = 0;
        else
            xu(i, j) = 1;
        end
    end
end
% STEP 4: Forming the deactivation vector for every time series
for i = 1:size(trs_n, 1)
    trs = trs_n(i, :)
    for j = 1: length(trs)
        if trs(1, j)>l
            xl(i, j) = 0;
        else
            xl(i, j) = -1;
        end
    end
end
% STEP 5: Calculation of accordance and discordance values
% Version 2
for i = 1:size(trs_n, 1)
    for j = 1:size(trs_n, 1)
        E(i, j) = (sqrt((xu(i, :)*xu(i, :)' + xl(i, :)*xl(i, :)'))*sqrt((xu(j, :)*xu(j, :)' + xl(j, :)*xl(j, :)')));
        acc(i, j) = (xu(i, :)*xu(j, :)' + xl(i, :)*xl(j, :)')/E(i, j);
        dcc(i, j) = (xu(i, :)*xl(j, :)' + xl(i, :)*xu(j, :)')/E(i, j);
    end
end

acc_tf = issymmetric(acc);
dcc_tf = issymmetric(dcc);
ans1 = sqrt((xu(1, :)*xu(1, :)' + xl(1, :)*xl(1, :)'));
ans2 = sqrt((xu(1, :)*xu(1, :)' + xl(1, :)*xl(1, :)'));
a = xu(1, :)*xu(1, :)' + xl(1, :)*xl(1, :)';
E = ans1*ans2;
E = sqrt((xu(1, :).*xu(1, :) + (xl(1, :).*xl(1, :)))).*sqrt((xu(1, :).*xu(1, :) + (xl(1, :).*xl(1, :))));
E_t = sqrt((xu(1, :)*xu(1, :)' + (xl(1, :)*xl(1, :)'))).*sqrt((xu(1, :)*xu(1, :)' + (xl(1, :).*xl(1, :))));
acc = (xu(1, :).*xu(1, :)') + (xl(1, :).*xl(1, :))/E;
dcc = (xu(1, :).*xl(2, :)) + (xl(1, :).*xu(2, :))/E;
% Testing the accordance and discordance function
% a(z, z) = 1
azz = (xu(1, :)*xu(1, :)' + xl(1, :)*xl(1, :)')/(sqrt((xu(1, :)*xu(1, :)' + xl(1, :)*xl(1, :)'))*sqrt((xu(1, :)*xu(1, :)' + xl(1, :)*xl(1, :)')));
% d(z, z) = 0
dzz = (xu(1, :)*xl(1, :)' + xl(1, :)*xu(1, :)')/(sqrt((xu(1, :)*xu(1, :)' + xl(1, :)*xl(1, :)'))*sqrt((xu(1, :)*xu(1, :)' + xl(1, :)*xl(1, :)')));
% z(z, -z) = 0

% d(z, -z) = -1

% -------------------------------------------------------------------------------- %
%                    GENERATION OF ALL THE FC ESTIMATES
% -------------------------------------------------------------------------------- %
% Getting the list of all the subjects
dirpath = dir(rootdir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
% Loading the MNI template
cd('E:\LRCBH\Data\MNI_segmented');
m = load_untouch_nii('GMWM_mask.nii');
mask = single(m.img);
% Iterating throught the different functional regions
cd(atlas_loc);
files = dir('*.nii');
region_list = extractfield(files,'name')';
% GENERATION OF THE REGION * TIME SERIES DATA
j=1;
for j = 1:length(subloc)
    suboi = subloc{j};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    filelist = dir('*.nii');
    timeseries_voxel = [];
    timeseries_region = [];
    fprintf('Generating data for subject %s...\n',suboi);
    for i = 1:length(region_list)
        cd(atlas_loc);
        region = load_untouch_nii(region_list{i,1});
        roi = single(region.img);
        [x, y, z] = size(roi);
        roi_temp = reshape(roi, [1,x*y*z]);
        [~, indices] = find(roi_temp);
        for k=1:length(filelist)
            cd(current);
            vol = load_untouch_nii(filelist(k).name);
            I = vol.img;
            I_m = I.*mask;
            I_R = I_m.*roi;
            vt = reshape(I_R, [1,x*y*z]);
            %[~,idx] = find(vt);
            vt = vt(:,indices);
            data(k,:) = vt;
            avg_data = mean(data,2);
        end
        timeseries_voxel = double(horzcat(timeseries_voxel,data));
        timeseries_region = double(horzcat(timeseries_region,avg_data));
        voxel_loc{i,1} = indices;
        clear vol I I_R vt data avg_data region roi x y z roi_temp indices;
    end
    save(fullfile(data_savedir, sprintf('ROI_%s.mat',suboi)),'timeseries_voxel','timeseries_region','voxel_loc');
    clear timeseries_region timeseries_voxel voxel_loc;
end
cd('')
timeseries_region = timeseries_region(2:end, :)';
timeseries_region = timeseries_region(1:835, :);
r = corrcoef(timeseries_region);
y = 1:835;
plot(timeseries_region(1, :), y);
[acc_est, dcc_est, u, l] = accordance_discordance_estimate(timeseries_region', 0.8);
% Applying mask and saving a image
% Making sure the mask and the functional images are in the same space
clear all;
cd('E:\LRCBH\Data\MNI_segmented');
m = load_untouch_nii('GMWM_mask.nii');
mask = single(m.img);
ex_sub_loc = 'E:\LRCBH\Data\COBRE-MNI\Individual_data\Useful\008'
cd(ex_sub_loc);
filelist = dir('*.nii');
savepath = 'E:\LRCBH\Data\COBRE-MNI\Individual_data\Masked';
save_prefix = 'masked';
for i = 1:length(filelist)
    cd(ex_sub_loc);
    vol = load_untouch_nii(filelist(i).name);
    func = single(vol.img);
    res = func.*mask;
    vol.img = res;
    Sname = [filelist(i).name(1:16), '_' ,save_prefix, '.nii'];
    save_untouch_nii(vol,[savepath,'\',Sname]);
end
% Extracting time series from the masked images
cd(savepath);
filelist = dir('*.nii');
timeseries_voxel = [];
timeseries_region = [];
atlas_loc = 'E:\LRCBH\Atlas\Willard_with_overlap_2mm\1.Merged';
cd(atlas_loc);
files = dir('*.nii');
region_list = extractfield(files,'name')';
for i = 1:length(region_list)
    cd(atlas_loc);
    region = load_untouch_nii(region_list{i,1});
    roi = single(region.img);
    [x, y, z] = size(roi);
    roi_temp = reshape(roi, [1,x*y*z]);
    [~, indices] = find(roi_temp);
    for k=1:length(filelist)
        cd(savepath);
        vol = load_untouch_nii(filelist(k).name);
        I = vol.img;
        I_R = I.*roi;
        vt = reshape(I_R, [1,x*y*z]);
        %[~,idx] = find(vt);
        vt = vt(:,indices);
        data(k,:) = vt;
        avg_data = mean(data,2);
    end
    timeseries_voxel = double(horzcat(timeseries_voxel,data));
    timeseries_region = double(horzcat(timeseries_region,avg_data));
    voxel_loc{i,1} = indices;
    clear vol I I_R vt data avg_data region roi x y z roi_temp indices;
end
% -------------------------------------------------------------------------------- %
%                   TIME SERIES EXTRACTION : INDIVIDUAL RSN NETWORKS
% -------------------------------------------------------------------------------- %
% All the directory locations
pat_dir = 'W:\LRCBH\COBRE data';
atlas_dir = 'W:\Atlas\Willard shirer atlas\fROIs_90\2.Individual';
% Getting the list of all the subjects
dirpath = dir(pat_dir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
% Getting the different RSN groups
cd(atlas_dir);
atlaspath = dir(atlas_dir);
rsndir = [atlaspath(:).isdir];
rsngrp = {atlaspath(rsndir).name}';
rsngrp(ismember(rsngrp,{'.','..'})) = [];
clearvars -except pat_dir atlas_dir subloc rsngrp;
% GENERATION OF THE REGION * TIME SERIES DATA
for j = 1:length(subloc)
    suboi = subloc{j};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    filelist = dir('*.nii');
    timeseries_voxel = [];
    timeseries_region = [];
    fprintf('Generating data for subject %s...\n',suboi);
    for i = 1:length(region_list)
        cd(atlas_loc);
        region = load_nii(region_list{i,1});
        roi = single(region.img);
        [x, y, z] = size(roi);
        roi_temp = reshape(roi, [1,x*y*z]);
        [~, indices] = find(roi_temp);
        for k=1:length(filelist)
            cd(current);
            vol = load_untouch_nii(filelist(k).name);
            I = vol.img;
            I_R = I.*roi;
            vt = reshape(I_R, [1,x*y*z]);
            %[~,idx] = find(vt);
            vt = vt(:,indices);
            data(k,:) = vt;
            avg_data = mean(data,2);
        end
        timeseries_voxel = double(horzcat(timeseries_voxel,data));
        timeseries_region = double(horzcat(timeseries_region,avg_data));
        voxel_loc{i,1} = indices;
        clear vol I I_R vt data avg_data region roi x y z roi_temp indices;
    end
    save(fullfile(data_savedir, sprintf('ROI_%s.mat',suboi)),'timeseries_voxel','timeseries_region','voxel_loc');
    clear timeseries_region timeseries_voxel voxel_loc;
end


% Checking the orientation of the nii image
cd('E:\LRCBH\Data\COBRE-MNI\Individual_data\Masked');
im = load_untouch_nii('s8_MNI_012_00001_masked.nii');
im1 = load_nii('s8_MNI_012_00001_masked.nii');
m = 
m_rr = rri_orient(m);
im_rr = rri_orient(im);
cd('E:\LRCBH\Atlas\Willard_with_overlap_2mm\1.Merged');
mask = load_nii('2mm_merged_LECN.nii');
































