% clc
% clear
%============================================================
% skull strip T1 images, prepare for ants coregistration;
%==========================================================

clc
clear
[~,~,RAW] = xlsread('Z:\XZ\Cobre_Data_Quality.xlsx','tmp1');
% [~,~,RAW] = xlsread('Y:\COBRE_SCANS\excel_file_XZ\COBRE_all_data_rerun_09122019.xlsx','Year_2_need');
subject_list = RAW(2:end,1);
bold_id = cell2mat(RAW(2:end,7));
T1_id = cell2mat(RAW(2:end,9));
off_session_id = RAW(2:end,12); 
on_session_id = RAW(2:end,11);
Time_point=RAW(2:end,4);
session_id = RAW(2:end,13);

for i = 1:length(bold_id)
    if isnan(bold_id(i))
        continue;
    end
    script_1_realignment(subject_list(i),bold_id(i),off_session_id{i},on_session_id{i});
    script_2_1_prepare_distortion_correction_t(subject_list{i},bold_id(i),off_session_id{i},on_session_id{i},Time_point{i});
    script_2_2_split_4D_into_3D(subject_list(i),bold_id(i),off_session_id{i},on_session_id{i});
    script_3_segmentation_T1(subject_list(i),T1_id(i),session_id{i});
end

function script_3_segmentation_T1(subject_list,T1_id,session_id)
if isnan(session_id)
    return;
end
if isnumeric(session_id)
    session_id = num2str(session_id);
end

All_subject_path = 'Y:\COBRE_SCANS';
% All_result_path = 'Y:\COBRE_SCANS\processed_COBRE_scan';
% subject_list = {'LEVESQUE_DAVID_T','STOTT_VAR_C','SCOTT_JEFFREY','PAGONE_VINCENT_C'};
N_sub = length(subject_list);
for subject_index = 1:N_sub
    subject_path = [All_subject_path,'\',subject_list{subject_index}]
    if ~exist(subject_path,'dir')
        continue;
    end
    t1_name = dir([subject_path,'\',session_id(1),'*',sprintf('%03d',T1_id(1)),'*MPRAGE*']);
    t1_path = [subject_path,'\',t1_name(3).name];
    t1_file = selectImageDir(t1_path,'output.nii');
    SPM12_one_sub_segmentation(t1_path,'.nii');
    t1_mask = WMGMCSF_maskGeneration(t1_path);
    skullstrip(t1_file{1},t1_mask,t1_path,'t1_brain.nii');  
end



% Modifying the code to include all the ROI groups at once
% All the directory locations
pat_dir = 'F:\LRCBH\Data\COBRE_FLUENCY\fMRI\Useful';
atlas_dir = 'F:\LRCBH\Atlas\Willard shirer atlas\fROIs_90\2.Individual\Work';
savedir = 'F:\LRCBH\SSM data\GMWM';
% Loading the mask
cd('F:\LRCBH\Data\MNI_segmented');
m = load_nii('GMWM_mask.nii');
mask = single(m.img);
% Getting the list of all the subjects
dirpath = dir(pat_dir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
% Getting the different RSN groups
items = dir(atlas_dir);
rsn_names = {items.name}';
rsn_names(ismember(rsn_names,{'.','..'})) = [];
for rsn = 1:length(rsn_names)
    fprintf('Generating data for RSN %s...\n', rsn_names{rsn, 1});
    clearvars -except pat_dir atlas_dir savedir subloc roi_list mask rsn_names rsn;
    % Generating the time series data
    for pat = 1:length(subloc)
        suboi = subloc{pat};
        current_pat = strcat(pat_dir,'\', suboi);
        cd(current_pat);
        filelist = dir('*.nii');
        %fprintf('Generating data for subject %s...\n',suboi);
        cd([atlas_dir, '/',rsn_names{rsn, 1}]);
        files = dir('*.nii');
        region_list = extractfield(files,'name')';
        timeseries_voxel = [];
        for i = 1:length(region_list)
            cd([atlas_dir, '/',rsn_names{rsn, 1}]);
            region = load_nii(region_list{i,1});
            roi = single(region.img);
            roi_m = roi.*mask;
            [x, y, z] = size(roi_m);
            roi_temp = reshape(roi_m, [1,x*y*z]);
            [~, indices] = find(roi_temp);
            voxel_loc{1,1} = indices;
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
            end
            timeseries_voxel = double(horzcat(timeseries_voxel,data));
            clear vol I I_R vt data avg_data region roi x y z roi_temp indices;
            timeseries_voxel = timeseries_voxel';
            data_savedir = strcat(savedir,'\',rsn_names{rsn, 1}, '\','ROI_',num2str(i));
            if ~exist(data_savedir, 'dir')
                mkdir(data_savedir);
            end
            save(fullfile(data_savedir, sprintf('%s_roi%d_%s.mat',rsn_names{rsn, 1},i, suboi)),'timeseries_voxel','voxel_loc');
            clear timeseries_voxel voxel_loc;
            timeseries_voxel = [];
        end
    end
end
