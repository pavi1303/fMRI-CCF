% -------------------------------------------------------------------------------- %
%                   TIME SERIES EXTRACTION : INDIVIDUAL RSN NETWORKS
% -------------------------------------------------------------------------------- %
% All the directory locations
pat_dir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data\Useful';
atlas_dir = 'E:\LRCBH\Atlas\Willard shirer atlas\fROIs_90\2.Individual';
savedir = 'E:\LRCBH\Data\Voxel timeseries data\GMWM';
% Loading the mask
cd('E:\LRCBH\Data\MNI_segmented');
%m = load_nii('c1MNI152_T1_2mm.nii');
m = load_nii('GMWM_mask.nii');
mask = single(m.img);
% Getting the list of all the subjects
dirpath = dir(pat_dir);
subdir = [dirpath(:).isdir]; 
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
subloc(84:106, :) = [];
% Getting the different RSN groups
cd(atlas_dir);
atlaspath = dir(atlas_dir);
rsndir = [atlaspath(:).isdir];
rsngrp = {atlaspath(rsndir).name}';
rsngrp(ismember(rsngrp,{'.','..'})) = [];
clearvars -except pat_dir atlas_dir savedir subloc rsngrp mask;
% Generating the time series data
for pat = 1:length(subloc)
    suboi = subloc{pat};
    current_pat = strcat(pat_dir,'\', suboi);
    cd(current_pat);
    filelist = dir('*.nii');
    fprintf('Generating data for subject %s...\n',suboi);
    for rsn = 1:length(rsngrp)
        rsnoi = rsngrp{rsn};
        %fprintf('%s ROI\n', rsnoi);
        current_rsn = strcat(atlas_dir, '\', rsnoi);
        cd(current_rsn);
        files = dir('*.nii');
        region_list = extractfield(files,'name')';
        timeseries = [];
        timeseries_voxel = [];
        for i = 1:length(region_list)
            cd(current_rsn);
            region = load_nii(region_list{i,1});
            roi = single(region.img);
            roi_m = roi.*mask;
            [x, y, z] = size(roi_m);
            roi_temp = reshape(roi_m, [1,x*y*z]);
            [~, indices] = find(roi_temp);
            voxel_loc{i,1} = indices;
            for k=1:length(filelist)
                cd(current_pat);
                vol = load_nii(filelist(k).name);
                I = vol.img;
                I_m = I.*mask;
                I_R = I_m.*roi_m;
                vt = reshape(I_R, [1,x*y*z]);
                %[~,idx] = find(vt);
                vt = vt(:,indices);
                data(k,:) = vt;
                avg_data = mean(data,2);
            end
            timeseries_voxel = double(horzcat(timeseries_voxel,data));
            timeseries = double(horzcat(timeseries,avg_data));
            clear vol I I_R vt data avg_data region roi x y z roi_temp indices;
        end
        timeseries = timeseries';
        timeserie_voxel = timeseries_voxel';
        data_savedir = strcat(savedir,'\',num2str(rsn), '.',rsnoi);
        if ~exist(data_savedir, 'dir')
            mkdir(data_savedir);
        end
        save(fullfile(data_savedir, sprintf('%s_roi_%s.mat',rsnoi, suboi)),'timeseries','timeseries_voxel','voxel_loc');
        clear timeseries voxel_loc;
    end
end
% -------------------------------------------------------------------------------- %
%               FUNCTIONAL CONNECTIVITY ESTIMATION : ALL TYPES
% -------------------------------------------------------------------------------- %
savedir = 'E:\LRCBH\Projects\COBRE\Home\Timeseries_data\GMWM';
fcn_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\Masked\GMWM\1.Functional connectivity';
cd(savedir);
atlaspath = dir(savedir);
rsndir = [atlaspath(:).isdir];
rsngrp = {atlaspath(rsndir).name}';
rsngrp(ismember(rsngrp,{'.','..'})) = [];
for rsn = 1:length(rsngrp)
    rsnoi = rsngrp{rsn};
    dirpath = dir(strcat(savedir, '\', rsnoi));
    subpath = {dirpath.name}';
    subpath(ismember(subpath,{'.','..'})) = [];
    for sub = 1:length(subpath)
        suboi = subpath{sub};
        suboi_temp = reverse(suboi);
        tc_data = load(strcat(savedir, '\', rsnoi,'\',suboi),'timeseries');
        tc = (tc_data.timeseries);
        % Full correlation
        fcn_full_pearson = corr(tc(:, 16:850)', 'type', 'Pearson'); % Full correlation - Pearson
        fcn_full_spearman = corr(tc(:, 16:850)', 'type','Spearman'); % Full correlation - Spearmann
        fcn_full_kendall = corr(tc(:, 16:850)', 'type', 'Kendall'); % Full correlation - Kendall Tau
        if ~exist(strcat(fcn_savedir,'\Full correlation - Pearson','\',rsnoi), 'dir')
            mkdir(strcat(fcn_savedir,'\Full correlation - Pearson','\',rsnoi));
        end 
        if ~exist(strcat(fcn_savedir,'\Full correlation - Spearman','\',rsnoi), 'dir')
            mkdir(strcat(fcn_savedir,'\Full correlation - Spearman','\',rsnoi));
        end
        if ~exist(strcat(fcn_savedir,'\Full correlation - Kendall','\',rsnoi), 'dir')
            mkdir(strcat(fcn_savedir,'\Full correlation - Kendall','\',rsnoi));
        end
        save(fullfile(strcat(fcn_savedir,'\Full correlation - Pearson','\',rsnoi), sprintf('FCN_full_pearson_%s',reverse(suboi_temp(5:7)))),'fcn_full_pearson');
        save(fullfile(strcat(fcn_savedir,'\Full correlation - Spearman','\',rsnoi), sprintf('FCN_full__spearman_%s',reverse(suboi_temp(5:7)))),'fcn_full_spearman');
        save(fullfile(strcat(fcn_savedir,'\Full correlation - Kendall','\',rsnoi), sprintf('FCN_full_kendall_%s',reverse(suboi_temp(5:7)))),'fcn_full_kendall');
        % Partial correlation
        fcn_partial = partialcor(tc(:, 16:850)'); % Partial correlation
        if ~exist(strcat(fcn_savedir,'\Partial correlation','\',rsnoi), 'dir')
            mkdir(strcat(fcn_savedir,'\Partial correlation','\',rsnoi));
        end
        save(fullfile(strcat(fcn_savedir,'\Partial correlation','\',rsnoi), sprintf('FCN_partial_%s',reverse(suboi_temp(5:7)))),'fcn_partial');
        % Accordance and discordance
        [fcn_acc, fcn_dcc, ~, ~] = accordance_discordance_estimate(tc(:, 16:850), 0.80);
        if ~exist(strcat(fcn_savedir,'\Accordance','\',rsnoi), 'dir')
            mkdir(strcat(fcn_savedir,'\Accordance','\',rsnoi));
        end
        if ~exist(strcat(fcn_savedir,'\Discordance','\',rsnoi), 'dir')
            mkdir(strcat(fcn_savedir,'\Discordance','\',rsnoi));
        end
        save(fullfile(strcat(fcn_savedir,'\Accordance','\',rsnoi), sprintf('FCN_accordance_%s',reverse(suboi_temp(5:7)))),'fcn_acc');
        save(fullfile(strcat(fcn_savedir,'\Discordance','\',rsnoi), sprintf('FCN_discordance_%s',reverse(suboi_temp(5:7)))),'fcn_dcc');
    end
end
% -------------------------------------------------------------------------------- %
%                       GENERATION OF THE TARGET VARIABLE
% -------------------------------------------------------------------------------- %
cd('E:\LRCBH\Projects\COBRE\Results\Documents\Excel');
order = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 13 107 13]);
fcn_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\Masked\GM\1.Functional connectivity';
regdata_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\Masked\GM\2.Regression data';
dirloc = dir(fcn_savedir);
fcntype = get_dir(fcn_savedir);
% Getting the individual rsn names
rsngrp = get_dir(strcat(fcn_savedir,'\',fcntype{1,1}));
rsngrp_temp = rsngrp;
for k = 1 : length(rsngrp_temp)
    cellContents = rsngrp_temp{k};
    rsngrp_temp{k} = cellContents(3:end);
end
for k = 2 : 6
    cellContents = rsngrp_temp{k};
    rsngrp_temp{k} = cellContents(2:end);
end
% Need to import the design matrix and include that when saving the
designmat_path  = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis\Data';
cd(designmat_path);
data = load("regdata_merged_MI_own.mat");
X = data.X(:, 1:6);
cd('E:\LRCBH\Projects\COBRE\Results\Documents\Excel');
X(:,4) = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 9 107 9]);
X(:,1) = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 10 107 10]);
X(:, 3) = X(:, 1).*X(:, 2);
%X_mci = data.X_mci;
%X_nc = data.X_nc;
% -------------------------------------------------------------------------------- %
%                           ACCORDANCE AND DISCORDANCE
% -------------------------------------------------------------------------------- %
fcntypeoi = fcntype{1,1};
rsngrp = get_dir(strcat(fcn_savedir,'\',fcntypeoi));
for j = 1: length(rsngrp)
    cd(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j}));
    dirloc = dir(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j}));
    subloc = {dirloc.name}';
    subloc(ismember(subloc,{'.','..'})) = [];
    subloc = subloc(order);
    for k = 1:length(subloc)
        sub_fcn_acc = load(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j},'\',subloc{k}),'fcn_acc');
        sub_fcn_acc = sub_fcn_acc.fcn_acc;
        sub_fcn_dcc = load(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j},'\',subloc{k}),'fcn_dcc');
        sub_fcn_dcc = sub_fcn_dcc.fcn_dcc;
        Y_acc(k, :) = getfcn_vec(sub_fcn_acc);
        Y_dcc(k, :) = getfcn_vec(sub_fcn_dcc);
    end
    if ~exist(strcat(regdata_savedir,'\',fcntypeoi,'\Accordance'))
        mkdir(strcat(regdata_savedir,'\',fcntypeoi,'\Accordance'))
    end
    if ~exist(strcat(regdata_savedir,'\',fcntypeoi,'\Discordance'))
        mkdir(strcat(regdata_savedir,'\',fcntypeoi,'\Discordance'))
    end
    save(fullfile(strcat(regdata_savedir,'\',fcntypeoi,'\Accordance'), sprintf('regdata_acc_%s.mat',rsngrp_temp{j})), 'Y_acc', 'X');
    save(fullfile(strcat(regdata_savedir,'\',fcntypeoi,'\Discordance'), sprintf('regdata__dcc_%s.mat',rsngrp_temp{j})), 'Y_dcc', 'X');
    clear Y_acc Y_dcc
end
% -------------------------------------------------------------------------------- %
%                                         FULL CORRELATION
% -------------------------------------------------------------------------------- %
fcntypeoi = fcntype{2,1};
rsngrp = get_dir(strcat(fcn_savedir,'\',fcntypeoi));
for j = 1: length(rsngrp)
    cd(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j}));
    dirloc = dir(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j}));
    subloc = {dirloc.name}';
    subloc(ismember(subloc,{'.','..'})) = [];
    subloc = subloc(order);
    for k = 1:length(subloc)
        sub_fcn_pearson = load(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j},'\',subloc{k}),'fcn_full_pearson');
        sub_fcn_pearson = sub_fcn_pearson.fcn_full_pearson;
        sub_fcn_spearman = load(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j},'\',subloc{k}),'fcn_full_spearman');
        sub_fcn_spearman = sub_fcn_spearman.fcn_full_spearman;
        sub_fcn_kendall = load(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j},'\',subloc{k}),'fcn_full_kendall');
        sub_fcn_kendall = sub_fcn_kendall.fcn_full_kendall;
        Y_full_pearson(k, :) = getfcn_vec(sub_fcn_pearson);
        Y_full_spearman(k, :) = getfcn_vec(sub_fcn_spearman);
        Y_full_kendall(k, :) = getfcn_vec(sub_fcn_kendall);
    end
    if ~exist(strcat(regdata_savedir,'\',fcntypeoi,'\Pearson'))
        mkdir(strcat(regdata_savedir,'\',fcntypeoi,'\Pearson'))
    end
    if ~exist(strcat(regdata_savedir,'\',fcntypeoi,'\Spearman'))
        mkdir(strcat(regdata_savedir,'\',fcntypeoi,'\Spearman'))
    end
    if ~exist(strcat(regdata_savedir,'\',fcntypeoi,'\Kendall'))
        mkdir(strcat(regdata_savedir,'\',fcntypeoi,'\Kendall'))
    end
    save(fullfile(strcat(regdata_savedir,'\',fcntypeoi,'\Pearson'), sprintf('regdata_pearson_%s.mat',rsngrp_temp{j})), 'Y_full_pearson', 'X');
    save(fullfile(strcat(regdata_savedir,'\',fcntypeoi,'\Spearman'), sprintf('regdata_spearman_%s.mat',rsngrp_temp{j})), 'Y_full_spearman', 'X');
    save(fullfile(strcat(regdata_savedir,'\',fcntypeoi,'\Kendall'), sprintf('regdata_kendall_%s.mat',rsngrp_temp{j})), 'Y_full_kendall',  'X');
    clear Y_full_kendall Y_full_pearson Y_full_spearman
end
% -------------------------------------------------------------------------------- %
%                                         PARTIAL CORRELATION
% -------------------------------------------------------------------------------- %
fcntypeoi = fcntype{3,1};
rsngrp = get_dir(strcat(fcn_savedir,'\',fcntypeoi));
for j = 1:length(rsngrp)
    cd(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j}));
    dirloc = dir(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j}));
    subloc = {dirloc.name}';
    subloc(ismember(subloc,{'.','..'})) = [];
    subloc = subloc(order);
    for k = 1:length(subloc)
        sub_fcn_partial = load(strcat(fcn_savedir,'\',fcntypeoi,'\',rsngrp{j},'\',subloc{k}),'fcn_partial');
        sub_fcn_partial = sub_fcn_partial.fcn_partial;
        Y_partial(k, :) = getfcn_vec(sub_fcn_partial);
    end
    if ~exist(strcat(regdata_savedir,'\',fcntypeoi))
        mkdir(strcat(regdata_savedir,'\',fcntypeoi))
    end
    save(fullfile(strcat(regdata_savedir,'\',fcntypeoi), sprintf('regdata_partial_%s.mat',rsngrp_temp{j})), 'Y_partial', 'X');
    clear Y_partial
end
% -------------------------------------------------------------------------------- %
%             MEAN FUNCTIONAL CONNECTIVITY MATRICES FORMATION
% -------------------------------------------------------------------------------- %
cd('E:\LRCBH\Projects\COBRE\Results\Documents\Excel');
order = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 13 107 13]);
fcn_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\Masked\GMWM\1.Functional connectivity';
fcntype = get_dir(fcn_savedir);
for i = 1:length(fcntype)
    rsngrp = get_dir(strcat(fcn_savedir, '\', fcntype{i}));
    for j = 1:length(rsngrp)
        cd(strcat(fcn_savedir, '\', fcntype{i}, '\', rsngrp{j}));
        dirloc = dir(strcat(fcn_savedir, '\', fcntype{i}, '\', rsngrp{j}));
        subloc = {dirloc.name}';
        subloc(ismember(subloc,{'.','..'})) = [];
        subloc = subloc(order);
        fc_nc = cell(1,51);
        fc_mci = cell(1, 55);
        for nc_idx = 1:length(fc_nc)
            sub_fcn = load(strcat(fcn_savedir, '\', fcntype{i}, '\', rsngrp{j},'\',subloc{nc_idx}));
            fc_nc(nc_idx) = struct2cell(sub_fcn);
        end
        nc_all = cat(3, fc_nc{:});
        fc_nc_mean = mean(nc_all, 3);
        clear sub_fcn;
        for mci_idx = 1:length(fc_mci)
            sub_fcn = load(strcat(fcn_savedir, '\', fcntype{i}, '\', rsngrp{j},'\',subloc{mci_idx}));
            fc_mci(mci_idx) = struct2cell(sub_fcn);
        end
        mci_all = cat(3, fc_mci{:});
        fc_mci_mean = mean(mci_all, 3);
        saveloc = strcat('E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\Masked\GMWM\2.Mean fc\', ...
            fcntype{i});
        if ~exist(saveloc)
            mkdir(saveloc)
        end
        save(fullfile(saveloc, sprintf('Mean_FC_%s',rsngrp{j,1}(4:end))),'fc_nc_mean', 'fc_mci_mean');
    end
end













