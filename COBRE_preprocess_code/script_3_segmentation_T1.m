% clc
% clear
%============================================================
% skull strip T1 images, prepare for ants coregistration;
%==========================================================
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
    t1_name = dir([subject_path,'\',session_id,'*',sprintf('%03d',T1_id),'*MPRAGE*']);
    t1_path = [subject_path,'\',t1_name(1).name];
    t1_file = selectImageDir(t1_path,'output.nii');
    SPM12_one_sub_segmentation(t1_path,'.nii');
    t1_mask = WMGMCSF_maskGeneration(t1_path);
    skullstrip(t1_file{1},t1_mask,t1_path,'t1_brain.nii');  
end




