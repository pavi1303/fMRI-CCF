% clc
% clear
%============================================================
% 1. split each AP PA se 4D image into 3D,
% 2. use the first 2 volume for distortion correction;
% 3. realgined the first 2 volume to mean EPI
% distortion correction file saved in processed_COBRE_scan
%==========================================================

function script_2_1_prepare_distortion_correction_t(sub_ID,bold_id,off_session_id,on_session_id,timepoint)
if sum(isnan(on_session_id))>0
    flag_on = 0;
    display('OFF run');
else
    flag_on = 1;
    display('ON run');
end
if isnumeric(on_session_id)
    on_session_id = num2str(on_session_id);
end
if isnumeric(off_session_id)
    off_session_id = num2str(off_session_id);
end
All_subject_path = 'Y:\COBRE_SCANS';
All_result_path = 'Y:\COBRE_SCANS\processed_COBRE_scan';
subject_path = [All_subject_path,'\',sub_ID]

subject_result_path = [All_result_path,'\',sub_ID,'\',timepoint,'\rsfMRI'];
if ~exist(subject_result_path,'dir')
    mkdir(subject_result_path);
end
if flag_on == 0
    rsfMRI_run_name =  ls([subject_path,'\',off_session_id,'*',sprintf('%03d',bold_id),'*bold*']);
    %                 se_PA_run_name = ls([subject_path,'\',off_session_id,'*se_AP*']);  mistake found 09122019
    %                 se_AP_run_name = ls([subject_path,'\',off_session_id,'*se_PA*']);
    se_PA_run_name = ls([subject_path,'\',off_session_id,'*se_PA*']);
    se_AP_run_name = ls([subject_path,'\',off_session_id,'*se_AP*']);
    run_number = [off_session_id,'_',sprintf('%03d',bold_id)];
else
    rsfMRI_run_name =  ls([subject_path,'\',on_session_id,'*',sprintf('%03d',bold_id),'*bold*']);
    %                 se_PA_run_name = ls([subject_path,'\',on_session_id,'*se_AP*']);
    %                 se_AP_run_name = ls([subject_path,'\',on_session_id,'*se_PA*']);
    se_PA_run_name = ls([subject_path,'\',on_session_id,'*se_PA*']);
    se_AP_run_name = ls([subject_path,'\',on_session_id,'*se_AP*']);
    run_number = [on_session_id,'_',sprintf('%03d',bold_id)];
end
distortion_correction_result_path = [subject_result_path,'\',run_number,'_distortion_correction'];
if ~exist(distortion_correction_result_path,'dir')
    mkdir(distortion_correction_result_path);
end
run_mean_file = selectImageDir([subject_path,'\',rsfMRI_run_name],'mean*.nii');

if ~isempty(se_PA_run_name) && ~isempty(se_AP_run_name)
    RAbatch = load('m_files\SPM_batch\realign_2_first.mat');
    nii_file_input_spm = cell(5,1);
    nii_file_input_spm{1} = run_mean_file{1};
    se_PA_file_t = selectImageDir([subject_path,'\',se_PA_run_name],'*.nii');
    for j = 2:3
        nii_file_input_spm{j,1} = [se_PA_file_t{1},',',num2str(j-1)];
    end
    se_AP_file_t = selectImageDir([subject_path,'\',se_AP_run_name],'*.nii');
    for j = 4:5
        nii_file_input_spm{j,1} = [se_AP_file_t{1},',',num2str(j-3)];
    end
    RAbatch.matlabbatch{1,1}.spm.spatial.realign.estwrite.data{1,1} = nii_file_input_spm;
    spm_jobman('run',RAbatch.matlabbatch);
    se_PA_file_t = selectImageDir([subject_path,'\',se_PA_run_name],'r*.nii');
    abk_4Dto3D(se_PA_file_t{1});
    se_AP_file_t = selectImageDir([subject_path,'\',se_AP_run_name],'r*.nii');
    abk_4Dto3D(se_AP_file_t{1});
    realgined_se_PA_file{1} = selectImageDir([subject_path,'\',se_PA_run_name],'r*001.nii');
    realgined_se_AP_file{1} = selectImageDir([subject_path,'\',se_AP_run_name],'r*001.nii');
    realgined_se_PA_file{2} = selectImageDir([subject_path,'\',se_PA_run_name],'r*002.nii');
    realgined_se_AP_file{2} = selectImageDir([subject_path,'\',se_AP_run_name],'r*002.nii');
end
for k = 1:2
    copyfile(realgined_se_AP_file{k}{1},[distortion_correction_result_path,'\r_se_AP_',sprintf('%03d',k),'.nii']);
    copyfile(realgined_se_PA_file{k}{1},[distortion_correction_result_path,'\r_se_PA_',sprintf('%03d',k),'.nii']);
end
