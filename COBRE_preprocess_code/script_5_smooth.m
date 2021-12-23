% clc
% clear
function script_5_smooth(subject_list,bold_id,session_id)
All_subject_path = 'Y:\COBRE_SCANS';
SM_fwhm = 8;
if isnumeric(session_id)
    session_id = num2str(session_id);
end
% spm_jobman('initcfg');
subject_path = [All_subject_path,'\',subject_list{1}];

%%%%%%%%%%%realignment
N_run = length(bold_id);
for ii = 1:N_run
    run_number = sprintf('%03d',bold_id(ii));
    rsfMRI_run_name = ls([subject_path,'\',session_id,'*',run_number,'*bold*']);
    run_path = [subject_path,'\',rsfMRI_run_name,'\intermediate'];
    if ~exist(run_path,'dir')
        continue;
    end
    nii_file = selectImageDir(run_path,'MNI_*.nii');
    SMbatch = load('m_files\SPM_batch\smooth.mat');
    SMbatch.matlabbatch{1,1}.spm.spatial.smooth.data = nii_file;
    SMbatch.matlabbatch{1,1}.spm.spatial.smooth.fwhm = [SM_fwhm,SM_fwhm,SM_fwhm];
    SMbatch.matlabbatch{1,1}.spm.spatial.smooth.prefix = ['s',num2str(SM_fwhm),'_'];
    spm_jobman('run',SMbatch.matlabbatch);
end
