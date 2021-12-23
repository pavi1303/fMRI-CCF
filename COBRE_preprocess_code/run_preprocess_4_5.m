clc
clear
% Time_point should be Baseline or Year_1; 

% subject_list = {'121','121','122','122','123'}; %first off then on;
% bold_id = [25,12,8,18,12];
% T1_id = [22,5,5,6,6];
% Time_point= {'Baseline','Baseline','Baseline','Baseline','Baseline','Baseline'};
% session_id ={'107382481','107366674','107396974','107396973','107361208'};


[NUM,TXT,RAW] = xlsread('Y:\XZ\Cobre_Data_Quality.xlsx','tmp1');
subject_list = RAW(2:end,1);
bold_id = NUM(:,1);
T1_id = NUM(:,3);
Time_point=RAW(2:end,4);
session_id = RAW(2:end,13);

% [NUM,TXT,RAW] = xlsread('Y:\COBRE_SCANS\excel_file_XZ\COBRE_all_data_rerun_09122019.xlsx','Year_2_need');
% subject_list = RAW(2:end,1);
% bold_id = NUM(:,3);
% T1_id = NUM(:,5);
% Time_point=RAW(2:end,4);
% session_id = RAW(2:end,13);

for i = 1:length(bold_id)
    if isnan(bold_id(i))
        continue;
    end
    if isnan(T1_id(i))
        continue;
    end
%     if i>3
        script_4_ants_coregistration_aMCI_fMRI_2_T1_2_MNI(subject_list(i),bold_id(i),T1_id(i),session_id{i},Time_point{i});
%     end
    script_5_smooth(subject_list(i),bold_id(i),session_id{i});
end