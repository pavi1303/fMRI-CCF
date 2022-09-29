%------------------------------------------------------------------------------%
%          WILLARD ROI ANALYSIS - FUNCTIONAL CONNECTIVITY
%------------------------------------------------------------------------------%
clc;
clear all;
close all;

% The timeseries and hence the FC metrics are ordered, do not forget to reorder the other independent
% variables during the multiple linear regression.
root_dir = 'F:\LRCBH\SSM data\GMWM\Language';
trs_savedir = 'F:\LRCBH\Willard_timeseries_data';
excel_dir = 'F:\LRCBH\Projects\COBRE\Results\Documents\Excel';
cd(excel_dir);
order = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 17 90 17]);
d = dir(root_dir);
isub = [d(:).isdir];
subdirs = {d(isub).name}';
subdirs(ismember(subdirs,{'.','..'})) = [];

for i = 1:numel(subdirs)
    cd(fullfile(root_dir, subdirs{i}));
    files = dir(fullfile(fullfile(root_dir, subdirs{i}), '*.mat'));
    filenames = {files.name}';
    mean_trs = zeros(size(filenames, 1), 850);
    for j = 1:length(filenames)
        sub_ts = load(strcat(fullfile(root_dir, subdirs{i}),'\',filenames{j}),'timeseries_voxel');
        mean_trs(j, :) = mean(sub_ts.timeseries_voxel);
    end
    mean_timeseries = mean_trs(order, :);
    saveloc = fullfile(trs_savedir, extractAfter(root_dir, 23));
    if ~exist(saveloc,'dir')
        mkdir(saveloc);
    end
    save(fullfile(saveloc, sprintf('Mean_timeseries_%s_%s.mat',extractAfter(root_dir, 23), subdirs{i})),'mean_timeseries');
    clear mean_timeseries mean_trs files filenames;
end

% Using the timeseries data and deriving the functional connectivity
% matrices
root_dir = 'F:\LRCBH\Timeseries_data_willard_ROI\Language';
d = dir(root_dir);
isub = [d(:).isdir];
subdirs = {d(isub).name}';
subdirs(ismember(subdirs,{'.','..'})) = [];











