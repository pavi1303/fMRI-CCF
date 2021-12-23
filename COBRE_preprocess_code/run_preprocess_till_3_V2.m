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