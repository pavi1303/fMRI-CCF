% clc
% clear
%=====================================================================
% split realigned 4D image into 3D images, prepare for distortion
% correction, because Applytopup on 4D not stable; saved in
% subject_path\intermediate;
%====================================================================
function script_2_2_split_4D_into_3D(subject_list,bold_id,off_session_id,on_session_id)
if isnan(on_session_id)
    flag_on = 0;
    display('OFF run');
else
    flag_on = 1;
    display('ON run');
end
if isnumeric(off_session_id)
    off_session_id = num2str(off_session_id);
end
if isnumeric(on_session_id)
    on_session_id = num2str(on_session_id);
end
All_subject_path = 'Y:\COBRE_SCANS';
% subject_list = {'LEVESQUE_DAVID_T','PAGONE_VINCENT_C','SCOTT_JEFFREY','STOTT_VAR_C'};
N_sub = length(subject_list);
for subject_index = 1:N_sub
    subject_path = [All_subject_path,'\',subject_list{subject_index}]
    if ~exist(subject_path,'dir')
        continue;
    end
    N_run = length(bold_id);
    if N_run ~= 0
        for i = 1:N_run
            run_number = sprintf('%03d',bold_id(i));
            if flag_on == 0
                rsfMRI_run_name = ls([subject_path,'\',off_session_id,'*',sprintf('%03d',bold_id(i)),'*bold*']);
            else
                rsfMRI_run_name = ls([subject_path,'\',on_session_id,'*',sprintf('%03d',bold_id(i)),'*bold*']);
            end
            run_path = [subject_path,'\',rsfMRI_run_name];
            save_path_3D = [run_path,'\intermediate'];
            if ~exist(save_path_3D,'dir')
                mkdir(save_path_3D);
            end
%             realgined_file = selectImageDir(run_path,['r*',run_number,'*.nii']);
            realgined_file = selectImageDir(run_path,['r*.nii']);
            spm_file_split(realgined_file{1}, save_path_3D);
            delete(realgined_file{1});
        end
    end
end