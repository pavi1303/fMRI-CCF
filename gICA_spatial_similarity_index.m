clc
close all;
clear all;
%------------------------------------------------------------------------------%
%                     IMPORTING LIST OF USEFUL DIRECTORIES
%------------------------------------------------------------------------------%

rootdir = 'W:\LRCBH\COBRE_rsfMRI\fMRI\Useful';
pca_savedir = 'W:\ICA_similarity\PCA';
ica_savedir = 'W:\ICA_similarity\ICA';
dr_savedir = 'W:\ICA_similarity\Dual regression';
trs_savedir = 'W:\ICA_similarity\Timeseries';
tcat_savedir = 'W:\ICA_similarity\Temporal concatenation';

dirpath = dir(rootdir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];

%------------------------------------------------------------------------------%
%                     ACQUIRING CRUCIAL VOXEL INFORMATION
%------------------------------------------------------------------------------%

cd('W:\LRCBH\MNI_segmented');
%cd('W:\MNI_segmented')
m = load_untouch_nii('standard_binary.nii');
M = m.img;
nii_temp = m;
[x, y, z] = size(M);
mask_temp = reshape(M, [1,x*y*z]);
[~, indices] = find(mask_temp);
trs = 850;
n_pca = 200;

%------------------------------------------------------------------------------%
%                                 PCA REDUCTION
%------------------------------------------------------------------------------%

% Iterating through each of the subjects
for sub=1:length(subloc)
    suboi = subloc{sub};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    files = dir('*.nii');
    vt_data = zeros(trs, 228453);
    %Generating the voxel*timeseries
    for j=1:length(files)
        i = load_untouch_nii(files(j).name);
        I = i.img;
        I_M = I.*M;
        vt = reshape(I_M, [1,x*y*z]);
        [~,idx] = find(vt);
        vt = vt(:,idx);
        vt_data(j,:) = vt;
    end
    cd(trs_savedir);
    save(fullfile(trs_savedir, sprintf('VT_%s.mat',suboi)),'vt_data','-v7.3');
    clear I I_M vt;
    %PCA reduction
    c = cov(vt_data');
    [eigvec,eigenval,~]=pcacov(c);
    EV=eigvec(:,1:n_pca);
    pca_data=EV'*vt_data;
    cd(pca_savedir);
    save(fullfile(pca_savedir, sprintf('PCA_%s.mat',suboi)),'pca_data','-v7.3');
    clear pca_data;
end


%------------------------------------------------------------------------------%
%         TEMPORAL CONCATENATION, SPATIAL ICA - TIMESERIES
%------------------------------------------------------------------------------%

tcat_data = temporal_concat(subloc,trs_savedir);
save(fullfile(tcat_savedir, sprintf('tcat_data.mat',suboi)),'tcat_data','-v7.3');

clearvars -except tcat_data ica_savedir

method=11;
a1=1;
var_normal=1;
eps=1E-6;
A0=[];
shift=[];
Sigma=0;
determine_flip=1;
npca=50;

[S,W,White,E,eigval,convergence,A,B,A_reduced,X_reduced,Sigma_reduced]=...
    ica_DC_improved(tcat_data,Sigma,method,eps,npca,A0,a1,var_normal,shift,determine_flip);
save(fullfile(ica_savedir,sprintf('gica_%d_result.mat',npca)),'S','W','White','A','-v7.3');

%------------------------------------------------------------------------------%
%                  GENERATING NII VERSION OF ICA COMPONENTS
%------------------------------------------------------------------------------%

cd(ica_savedir);
ica_dat = load('gICA_50_result.mat','S');
S = ica_dat.S;
S = double(S);
save_ica_nii(S,x,y,z,indices,m,'gICA',ica_savedir);

%------------------------------------------------------------------------------%
%                                        DUAL REGRESSION
%------------------------------------------------------------------------------%

dirloc = dir(trs_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(trs_savedir,'\', suboi);
    sub_data = load(current,'vt_data');
    sub_data = (sub_data.vt_data)';
    %Performing spatial regression
    ss_tc = (pinv(S))'*sub_data;
    %Performing temporal regression
    ss_sm = sub_data*pinv(ss_tc);
    ss_dr_dir = strcat(dr_savedir,'\',suboi(5:7));
    if ~exist(ss_dr_dir,'dir')
        mkdir(ss_dr_dir);
    end
    save(fullfile(ss_dr_dir,sprintf('dualregression.mat')),'ss_tc',"ss_sm");
    % Saving the independent components as nii images
    %save_ica_nii(ss_sm',x,y,z,indices,m,'ica_',ss_dr_dir);
end
fprintf('Dual regression done.\n')

%------------------------------------------------------------------------------%
%                         SPATIAL SIMILARITY CALCULATION
%------------------------------------------------------------------------------%
cd(ica_savedir);
ica_dat =  load(['gICA_' num2str(npca) '_result.mat'], 'S');
S = ica_dat.S;
spa_sm = zeros(length(subloc),size(S,1));

for icacomp = 1:size(S, 1)
    for sub = 1:length(subloc)
        cd(strcat(dr_savedir,'\', subloc{sub}));
        dr = load('dualregression.mat','ss_sm');
        s_sm = dr.ss_sm;
        spa_sm(sub, icacomp) = dot(s_sm(icacomp, :), S(icacomp, :))/dot(S(icacomp, :), S(icacomp, :));
        clear dr s_sm;
    end
end






























