% clc
% clear
%% now off and on have different T1s
function script_4_ants_coregistration_aMCI_fMRI_2_T1_2_MNI(sub_ID,bold_id,T1_id,session_id,time_point)
All_subject_path = 'Y:/COBRE_SCANS';
All_result_path = 'Y:/COBRE_SCANS/processed_COBRE_scan';
MNIfile='Y:/XZ/templates/MNI152_2mm/standard.nii';
subject_path = [All_subject_path,'/',sub_ID{1}]
subject_result_path = [All_result_path,'/',sub_ID{1},'/',time_point,'/rsfMRI'];
if ~exist(subject_result_path,'dir')
    mkdir(subject_result_path);
end
if isnumeric(session_id)
    session_id = num2str(session_id);
end
ants_path = [subject_result_path,'/',session_id,'_ANTS_2MNI'];
if ~exist(ants_path,'dir')
    mkdir(ants_path);
end
T1_name = dir([subject_path,'/',session_id,'*_',sprintf('%03d',T1_id),'_*']);
% for i = 1:length(T1_name)
%     t1_path =  [subject_path,'/',T1_name(i).name];
%     t1_file = selectImageDir_ants_in_windows(t1_path,'t1_brain.nii');
%     if ~isempty(t1_file)
%         break;
%     else
%         SPM12_one_sub_segmentation(t1_path,'output');
%         WMGMCSF_maskGeneration(t1_path,'output','t1_brain');
%         t1_file = selectImageDir_ants_in_windows(t1_path,'t1_brain.nii');
% %         break;
%     end
% end

 t1_path =  [subject_path,'/',T1_name.name];
 t1_file = selectImageDir_ants_in_windows(t1_path,'t1_brain.nii');
 if isempty(t1_file)
     SPM12_one_sub_segmentation(t1_path,'output');
     WMGMCSF_maskGeneration(t1_path,'output','t1_brain');
     t1_file = selectImageDir_ants_in_windows(t1_path,'*t1_brain.nii');
 end
if ~exist([ants_path,'/SyNT12MNI1Warp.nii.gz'],'file')
    antsT12MNI=['antsRegistrationSyN.sh -d 3 -f ',MNIfile,' -m ',t1_file{1},' -o ',ants_path,'/SyNT12MNI'];
    system(['"C:/cygwin64/bin/bash" -c "',antsT12MNI,'"']);
end
N_run = length(bold_id);
if N_run ~= 0
    for i = 1:N_run
        run_number = sprintf('%03d',bold_id(i));
        run_name = ls([subject_path,'/',session_id,'*',run_number,'*bold*']);
        run_path = [subject_path,'/',run_name];
        run_mean_path = [run_path,'/corr_mean_file'];
        if exist(run_mean_path,'dir')
            rmdir(run_mean_path,'s');
        end
        mkdir(run_mean_path);
        %             mean_file = selectImageDir_ants_in_windows(run_mean_path,'corr_mean*.nii*');
        %             if isempty(mean_file)
        mean_file_before_corr = selectImageDir_ants_in_windows(run_path,'corr_*.nii.gz');
        if isempty(mean_file_before_corr)
            pause(3600);
            mean_file_before_corr = selectImageDir_ants_in_windows(run_path,'corr_*.nii.gz');
        end
        copyfile(mean_file_before_corr{1},run_mean_path);
        mean_file_before_corr = selectImageDir_ants_in_windows(run_mean_path,'corr_*.nii.gz');
        gunzip(mean_file_before_corr{1});
        delete(mean_file_before_corr{1});
        SPM12_one_sub_segmentation(run_mean_path,'corr_');
        WMGMCSF_maskGeneration(run_mean_path,'corr_','corr_mean');
        mean_file = selectImageDir_ants_in_windows(run_mean_path,'ST*corr_mean*.nii*');
        %             end
        antsfmri2T1=['ANTS 3 -m MI[',t1_file{1},',',mean_file{1},',1,32] -i 0 -o ',ants_path,'/',run_number,'_mean_2T1'];
        system(['"C:/cygwin64/bin/bash" -c "',antsfmri2T1,'"']);
        antsWarp=['WarpImageMultiTransform 3 ',mean_file{1},' ',ants_path,'/',run_number,'_mean_2MNIed.nii.gz -R ',MNIfile,...
            ' ',ants_path,'/SyNT12MNI1Warp.nii.gz ',ants_path,'/SyNT12MNI0GenericAffine.mat ',ants_path,'/',run_number,'_mean_2T1Affine.txt'];
        system(['"C:/cygwin64/bin/bash" -c "',antsWarp,'"']);
        run_all_bold_path = [run_path,'/intermediate'];
        run_file = selectImageDir_ants_in_windows(run_all_bold_path,'corr_r*');
        %             run_file_name = dir([run_all_bold_path,'/corr_r*']);
        parfor kk = 1:length(run_file)
            vol_id = sprintf('%05d',kk);
            antsWarp=['WarpImageMultiTransform 3 ',run_file{kk},' ',run_all_bold_path,'/MNI_',run_number,'_',vol_id,'.nii -R ',MNIfile,...
                ' ',ants_path,'/SyNT12MNI1Warp.nii.gz ',ants_path,'/SyNT12MNI0GenericAffine.mat ',ants_path,'/',run_number,'_mean_2T1Affine.txt'];
            system(['"C:/cygwin64/bin/bash" -c "',antsWarp,'"']);
        end
    end
end





