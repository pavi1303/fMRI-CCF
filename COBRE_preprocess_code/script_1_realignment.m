% % clc
% % clear
%============================================================
% realign 4D image to its mean;
% store .txt into motion file folder;
%==========================================================

function script_1_realignment(subject_list,bold_id,off_session_id,on_session_id)
if isnan(on_session_id)
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
NOV = 850;
for i = 1:length(subject_list)
    subject_path = [All_subject_path,'\',subject_list{i}];
    if ~exist(subject_path,'dir')
        continue;
    end
    %%%%%%%%%%%realignment
    N_run = length(bold_id);
    if N_run ~= 0
        for ii = 1:N_run
            if flag_on == 0
                rsfMRI_run_name = ls([subject_path,'\',off_session_id,'*',sprintf('%03d',bold_id(ii)),'*bold*']);
            else
                rsfMRI_run_name = ls([subject_path,'\',on_session_id,'*',sprintf('%03d',bold_id(ii)),'*bold*']);
            end
            RAbatch = load('m_files\SPM_batch\realign_2_mean.mat');
            run_path = [subject_path,'\',rsfMRI_run_name];
            nii_file = selectImageDir(run_path,'*.nii');
            nii_file_input_spm = cell(NOV,1);
            for j = 1:NOV
                nii_file_input_spm{j,1} = [nii_file{1},',',num2str(j)];
            end
            RAbatch.matlabbatch{1,1}.spm.spatial.realign.estwrite.data{1,1} = nii_file_input_spm;
            spm_jobman('run',RAbatch.matlabbatch);
            motion_path = [run_path,'\motion_file'];
            if ~exist(motion_path,'dir')
                mkdir(motion_path);
            end
            motion_file = selectImageDir(run_path,'r*.txt');
            for k = 1:length(motion_file)
                movefile(motion_file{k},motion_path);
            end
        end
    end
end