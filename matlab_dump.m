%Getting the required subdir
rootdir = 'E:\LRCBH\Data\COBRE-MNI\Trial';
savedir = 'E:\LRCBH\Results\Matlab\1.PCA';
dirpath = dir(rootdir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
%Loading the mask file
cd('E:\LRCBH\MNI_segmented');
m = load_untouch_nii('standard_binary.nii');
M = m.img;
[x, y, z] = size(M);
trs = 850;
%Iterating through each of the subjects
for i=1:length(subloc)
    suboi = subloc{i};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    files = dir('*.nii');
    temp = zeros(trs, 228453);
    for j=1:length(files)
        i = load_untouch_nii(files(j).name);
        I = i.img;
        I_M = I.*M;
        vt = reshape(I_M, [1,x*y*z]);
        [~,idx] = find(vt);
        vt = vt(:,idx);
        temp(j,:) = vt;
        clear vt, i, I;
    end
    vt_data = zscore(temp(16:end, :));
    clear temp;
    %Performing subject wise PCA reduction
    vt_data= vt_data - mean(vt_data,2);
    c = cov(vt_data');
    [V,D,explained]=pcacov(c);
    EV=V(:,1:200);
    PCA_red=EV'*vt_data;
    cd(savedir);
    save(fullfile(savedir, sprintf('PCA_%s.mat',suboi)),'PCA_red');
end

i=1;

files = dir('*.nii');
temp = zeros(trs, 228453);
for j = 1:length(files)
    i = load_untouch_nii(files(j).name);
    I = i.img;
    I_M = I.*M;
    vt = reshape(I_M, [1,x*y*z]);
    [~,idx] = find(vt);
    vt = vt(:,idx);
    temp(j,:) = vt;
    clear vt, i, I;
end
vt_data = zscore(temp(16:end, :));
clear temp
%Performing subject wise PCA reduction
vt_data= vt_data - mean(vt_data,2);
c = cov(vt_data');
[V,D,explained]=pcacov(c);
EV=V(:,1:200);
PCA_red=EV'*vt_data;
cd(savedir);
savename = strcat()
save(strcat(savedir,'\','PCA%s.mat',suboi),'PCA_red)';

savename = fullfile(savedir, sprintf('PCA%s.mat',suboi));
savename = fullfile(savedir, sprintf('PCA_%s.mat',suboi));
save(fullfile(savedir, sprintf('PCA_%s.mat',suboi)),'PCA_red');

baseFileName = sprintf('%s.mat', sixCenters{k});
fullMatFileName = fullfile(folder, baseFileName);
save(fullMatFileName, 'grouped');


root_dir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data';
subdir = [root_dir(:).isdir];
nameFolds = {root_dir(subdir).name}';
subfolder = string(nameFolds);
nameFolds(ismember(nameFolds,{'.','..'})) = [];
files = dir(fullfile(strcat(root_dir,'\',subfolder(1,1))), '*.nii');

addpath(genpath('E:\LRCBH\Data\COBRE-MNI\Trial\008'))
addpath(genpath('E:\LRCBH\MNI_segmented'))

data = load_untouch_nii('s8_MNI_012_00001.nii');
img_data = data.img;
rootdir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data';
filelist = dir(fullfile(rootdir, '**\*.*'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);

cd('C:\Users\PATTIAP\Desktop\Dataset')
data = load_untouch_nii('MNI-008_unmasked.nii');
mask = load_untouch_nii('standard_binary.nii');
mask_data = mask.img;
vt_data = xyz_to_q_ALT(data.img);
%Applying mask
masked_data = bsxfun(@times, img, mask_data);
[x,y,z] = size(mask_data);
masked_img = permute(masked_data, [4, 1, 2, 3]);
vt_data = reshape(masked_data, [1, x*y*z]);
vt_data_modified = nonzeros(vt_data)';
vt_data_modified = zscore(vt_data_modified);
ig_size = masked_img.
%%% For iterating through different folders
basePath = pwd;  %your base path which is in your case myTraining
allPaths = dir(basePath);  %get all directory content
subFolders = [allPaths(:).isdir]; %get only indices of folders
foldersNames = {allPaths(subFolders).name}'; % filter folders names
foldersNames(ismember(foldersNames,{'.','..'})) = []; %delete default paths for parents return '.','..'
for i=1:length(foldersNames), %loop through all folders
    tmp = foldersNames{i};  %get folder by index
    p = strcat([basePath '\']);
    currentPath =strcat([p tmp]); % add base to current folder
    cd(currentPath);   % change directory to new path
    files = dir('*.jpg'); % list all images in your path which in your case could be John or Mary
    for j=1:length(files), % loop through your images
        img = imread(files(j).name); % read each image and do what you want
    end
end

% My try 1 for masking the image from a patient subfolder
% Loading of the mask file way outside the for loop and getting its
% dimensions
cd('C:\Users\pavig\Downloads');
m = load_untouch_nii('standard_binary.nii');
M = m.img;
[xres, yres, zres] = size(M);
% Get the number of time series from the user
trs = 850;
% This is the inner for loop where the generation of the voxel time data
% happens and the conversion of the PCA also happens for each subject and
% then I save th result of PCA as a mat variable.
cd('W:\008');
files = dir('*.nii');
temp = zeros(trs, xres*yres*zres);
for j = 1:length(files)
    i = load_untouch_nii(files(j).name);
    I = i.img;
    I_M = I.*M;
    vt = reshape(I_M, [1,xres*yres*zres]);
    %vt = nonzeros(vt)';
    [~,idx] = find(vt);
    temp(j,:) = vt;
    clear vt, i, I;
end
idx = unique(v);
[~,v] = find(temp);
vt_data = temp(16:end, :);
X = vt_data;
%PCA analysis in matlab
X = X - mean(X,2);
c = cov(X');
[E,eigval,explained]=pcacov(c);
EE=E(:,1:200);
X_tilde=EE'*X;

mean_x=mean(X,2);
[coeff,score] = pca(vt_data);
reduced = score(



% You load the mask image outside the for loop
j =1;
data = load_untouch_nii(files(j).name);
img = data.img;
img_masked = img.*mask_data;
vt_data = reshape(I_M, [1, xres*yres*zres]);
vt = nonzeros(vt_data)';
vt = zscore(vt);

% Items to do
%[1] Load the ydata file from xz and make sure that the PCA process is being done
%correctly and also check the mask computation process and the indexes for
%it.
%[2] Include to save the voxeltime data in the file and then see if the
%indexes change - they shouldn't coz they are from the mask whose 0 and 1 locations are the same.
%[3] Use these indexes to then reconstruct the 4D array from the 2D data.
%[4] After you figure out the reconstruction part, then  prepare the code
%getting the individual slices of the 4D ICA data.
% [5] Test this out on the trial dataset. If it works - then now run the
% program and go to sleep. Check the results tomorrow.

fcn_savedir = 'E:\LRCBH\Results\Matlab\v2\4.FCN\All';
dirloc = dir(fcn_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(fcn_savedir,'\', suboi);
    sub_data = load(current,'fcn');
    sub_data = (sub_data.fcn)';
    corr_mat = tril(sub_data,-1);
    corr_vec = nonzeros(corr_mat);
    fnc_val{i,1} = corr_vec';
end
% Forming the Y data
Y = vertcat(fnc_val{:});
% Performing regression to obtain coefficients
for k=1:size(Y,2)
    lr_model{k,1} = fitlm(X(:,1:4),Y(:,k));
    Y_fitted(:,k) = lr_model{k,1}.Fitted;
end

for k = 1:size(Y,2)
    [b(:,k),~,~,~,stats(:,k)] = regress(Y(:,k),X_grp1);
end
p-val
[b1,~,~,~,stats1] = regress(Y(:,1),X_grp1);
mdl = fitlm(X_grp1,Y(:,1));
stats(3,:) = stats(3,:)/size(Y,2);
p = 0.05/size(Y,2);
idx = find(stats(3,:)<p);

%
[row,col] = find(corr_mat);
idx1 = horzcat(row,col);
sig_idx = idx1(idx,:);
% Inputs for the confound function
% Directory of the functional connectivity matrices
% The design matrix
[coeff,pval_corr,stats_corrected,n_compar,sig_loc]
S_grp1 =struct;
S_grp2 = struct;

%Using the regress function
[S_grp1.coeff, S_grp1.pval, S_grp1.stats, S_grp1.comparisons, S_grp1.sig_asso] = confound_sig(fcn_savedir,X_grp1,'fcn',corr_mat);
[S_grp2.coeff, S_grp2.pval, S_grp2.stats, S_grp2.comparisons, S_grp2.sig_asso] = confound_sig(fcn_savedir,X_grp2,'fcn',corr_mat);

% Generating the scatter plot (VF vs FC)
VF = repmat(fluency_ratio,[1 16]);
Y = grp_res.Yfitted_sig;
g = horzcat(repmat(1,[1,44]),repmat(0,[1,44]))';
h = gscatter(VF,Y,g);
gscatter(fluency_ratio,Y(:,1),g)

%
[f1,xi1] = ksdensity(pf);
[f2,xi2] = ksdensity(sf);
[f3,xi3] = ksdensity(fluency_ratio);

figure;
plot(xi1,f1)
lgd = legend('Phonemic fluency');
title('pd estimate - Phonemic fluency');
figure;
plot(xi2,f2)
lgd = legend('Semantic fluency');
title('pd estimate - Semantic fluency');
figure;
plot(xi3,f3)
lgd = legend('Fluency ratio');
title('pd estimate - Fluency ratio');

pval = 0.05;
pval_adjusted = 0.05/171;
dr_dir='E:\LRCBH\Results\Matlab\3.DR';
var_name = 'dualregression';
% Testing the assumptions for the multiple linear regression model
X = horzcat(regressor,interaction,covariates);
k=1;
lr_model{k,1} = fitlm(X,Y(:,k),'RobustOpts','ols');
interaction = [];
pval = rand(1,100);
sig_idx = rand(1,100);
details = ["voxel index","p-value"];
sig_voxel = vertcat(sig_idx,pval);
% Getting the Y variable
dr_dir = 'E:\LRCBH\Results\Matlab\3.DR\Unbiased'
var_name = 'dualregression';
rsn_no = 14;
dirloc = dir(dr_dir);
subdir = [dirloc(:).isdir];
subloc = {dirloc(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(dr_dir,'\', suboi);
    cd(current);
    sub_sm = load(var_name);
    sub_sm = (sub_sm.ss_sm)';
    sub_smoi(i,:) = sub_sm(rsn_no,:);
end
% Forming the Y data
Y = sub_smoi;
% Fitting the multiple linear regression model
% With the interaction term
[regress_withinter.X, regress_withinter.Y, regress_withinter.Yfitted, regress_withinter.residuals, regress_withinter.betas, regress_withinter.pvalue, ...
    regress_withinter.tstatistic, regress_withinter.alpha_level, regress_withinter.pval_interactionvar, ...
    regress_withinter.meanRsquared_original, regress_withinter.meanRsquared_adjusted, regress_withinter.sig_voxel] = ...
    regress_model('E:\LRCBH\Results\Matlab\3.DR\Unbiased',regressor,interaction, covariates,'dualregression',14,'with_interaction','E:\LRCBH\Results\Matlab\v2\5.Association');
% Without the interaction term
[regress_withinter.X, regress_withinter.Y, regress_withinter.fit_model, regress_withinter.betas, regress_withinter.pvalue, ...
    regress_withinter.tstatistic, regress_withinter.alpha_level, regress_withinter.pval_interactionvar, ...
    regress_withinter.meanRsquared_original, regress_withinter.meanRsquared_adjusted, regress_withinter.sig_voxel] = ...
    regress_model('E:\LRCBH\Results\Matlab\3.DR\Unbiased',regressor,[], covariates,'dualregression',14,'no_interaction','E:\LRCBH\Results\Matlab\v2\5.Association');


dr_dir='E:\LRCBH\Results\Matlab\3.DR\Unbiased';
var_name = 'dualregression';
rsn_no = 4;
dirloc = dir(dr_dir);
subdir = [dirloc(:).isdir];
subloc = {dirloc(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(dr_dir,'\', suboi);
    cd(current);
    sub_sm = load(var_name);
    sub_sm = (sub_sm.ss_sm)';
    sub_smoi(i,:) = sub_sm(rsn_no,:);
end

tstat_inter = tstat(:,4)';
tstat_sig = tstat_inter(:,sig_voxel(1,:));
sig_voxel(3,:) = tstat_sig;

% Generating the mean t-statistic maps for the significant voxels
j=1;
idx_grp = vox_grp{14,1}';
Y_sig_grp1 = mean(regress_result.Yfitted(1:44,idx_grp));
Y_sig_grp2 = mean(regress_result.Yfitted(45:end,idx_grp));

[~,~,~,stats] = ttest2(regress_result.Y(1:44,idx_grp),regress_result.Y(45:end,idx_grp));
temp = zeros(1,228453);
temp(1,idx_grp) = stats.tstat;
temp(1,idx_grp) = Y_sig_grp1;
save_ica_nii(temp,x,y,z,indices,m,'tstat_map','E:\LRCBH\Results\Matlab\v2\5.Association\Spatial_maps');

t_stat = tstat(:,2);
idx = find(t_stat>2.5);
details = horzcat(t_stat(idx,1),pval(idx,4));
details_mean = mean(details);



%----Finding the mean of the FCN matrix
path = 'E:\LRCBH\Results\Matlab\ICA_100_results\3.FCN';
dirloc = dir(path);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];



%%% Generating the overall and mean functional connectivity across the two groups
% OVERALL FUNCTIONAL CONNECTIVITY MATRIX
dirloc = dir(fcn_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
conn_thresh_pos = 0.5;
conn_thresh_neg = -0.5;
for i  = 1:length(subloc)
    sub = subloc{i};
    current = strcat(fcn_savedir,'\',sub);
    sub_data = load(current,'fcn');
    fcn_mat{1,i} = (sub_data.fcn);
end
X = cat(3,fcn_mat{:});
fcn_overall = mean(X,3);
L_mat = tril(fcn_overall,-1);
[row, col, ~] = find(L_mat>conn_thresh_pos | L_mat<conn_thresh_neg);
idx = horzcat(row, col);
val1 = zeros(50,50);
for i = 1:size(idx,1)
    val1(idx(i,1),idx(i,2)) = L_mat(idx(i,1),idx(i,2));
end
clearvars -except fcn_savedir, val1, idx, L_mat;
fcn(noise_idx,:) = 0;
fcn(:,noise_idx) = 0;
uncommon = setdiff(grp1.thresh_idx, grp2.thresh_idx);



f = rand(1,100)';
ratio = 1.677;
ratio_N = repmat(ratio, 100);
r = corrcoef(f, ratio);



y1 = 2 * randn(1,500);
y2 = 3 * randn(1,500) + 5;
y3 = 5 * randn(1,500) + 5;
y = [y1 y2 y3];
x = [ones(1,500) 2*ones(1,500) 3*ones(1,500)];
%------------------------------------------------------------------------------%
%                          CORRELATION WITH FLUENCY SCORE
%------------------------------------------------------------------------------%
%Based on function
%noise_idx = [7,10,11,12,23,24,26,33,37,34,39,40,42,45]; % ICA 50
noise_idx =[2,12,17,25,28,31,32,35,44,47,51,52,55,62,66,67,71,78,79,84,85,93,96,99]; % ICA 100
all_idx = [1:102]';
grp1_idx = [1,2,5,9,12,17,20,22,30,31,32,36,38,39,42,43,44,45,46,48,49,50,51,54,56,57,58,...
    59,60,61,62,64,67,69,72,75,76,78,80,81,82,83,84,86,87,88,89,94,100,101,102];
grp2_idx = setdiff(all_idx,grp1_idx)';
[all.subs,all.sub_cat,all.fcn_mean, all.L_mat1, all.thresh_idx, all.thresh_mat] = fcnmat_results(fcn_savedir, 0.5, -0.5, 50, all_idx, noise_idx);
[grp1.subs,grp1.sub_cat,grp1.fcn_mean, grp1.L_mat1, grp1.thresh_idx, grp1.thresh_mat] = fcnmat_results(fcn_savedir, 0.5, -0.5, 50, grp1_idx, noise_idx);
[grp2.subs,grp2.sub_cat,grp2.fcn_mean, grp2.L_mat1, grp2.thresh_idx, grp2.thresh_mat] = fcnmat_results(fcn_savedir, 0.5, -0.5, 50, grp2_idx, noise_idx);
save(fullfile('E:\LRCBH\Results\Matlab\ICA_100_results',sprintf('FCN_%d_results.mat',100)),'all','grp1','grp2','grp1_idx','grp2_idx');


%%%%% Finding the corresponding connection FC value in the other group
pos = find(grp2.rsn_loc(:,1)==70 & grp2.rsn_loc(:,2)==16);
cor = grp2.correlation_coefficient(pos,1);
%%% Finding the rsn location after restricting the analysis to the IFG, ITG
%%%etc
useful_rsn = [6,14,16,46,54,56,57,59,64,65,73,77,80,83,86,94];
[m,n] = ndgrid(useful_rsn,useful_rsn);
Z = [m(:),n(:)];
[~, Locb] = ismember(Z, fliplr(Z), 'rows');
idx = sort(Locb);
Z1 = flip(Z);
loc = find(Z(:,1)==Z(:,2));
Z(loc,:)=[];
loc1 = find(Z(:,:) == flip(Z(:,:)));
edge = (1:15:240)';
edge = (1:14:224)';
edge = (1:13:208)';
Z1 = Z(1:16,:);
Z2 = Z(edge,:);
Z(edge,:)=[];
pos = find(grp1.rsn_loc(:,1)==Z(1,1)& grp1.rsn_loc(:,2) ==Z(1,2));
for i = 1:size(Z,1)
    pos(i,:) = find(rsn_all(:,1)==Z(i,1) & rsn_all(:,2) ==Z(i,2));
end

pos1 = sort(pos);
corr_shortlisted = find(pos1<2850);
corr_shortlist = pos1(1:120,1);
rsn_oi = grp1.rsn_loc(corr_shortlist,:);
corr_comparison(1,:) = grp1.correlation_coefficient(corr_shortlist,1);
corr_comparison(2,:) = grp2.correlation_coefficient(corr_shortlist,1);
diff = minus(corr_comparison(1,:),corr_comparison(2,:));
diff_loc = find(0.3<diff | diff<-0.3);
diff_shortlisted = diff(1,diff_loc);
corr_shortlisted = corr_comparison(:,diff_loc);
rsn_shortlisted = rsn_oi(diff_loc,:);
%%%%
rsn = grp1.rsn_loc;
rsn1 = fliplr(rsn);
rsn_all = vertcat(rsn, rsn1);


%%%%%
cd('E:\LRCBH\Data\COBRE-MNI\Individual_data\Useful\008');
dat = load_untouch_nii('s8_MNI_012_00001.nii');
I=dat.img;
m = load_untouch_nii('PP264_all_ROIs_combined_1mm.nii');
M = m.img;
B = imresize3(M,[91 109 91]);% Power atlas combined 2 mm
[x, y, z] = size(B);
vt = reshape(I_M, [1,x*y*z]);
vt(isnan(vt))=0;
[~,idx] = find(vt);
vt = vt(:,idx);
%%%% Trying out the atlas
I_M = I.*B;
for j=1:length(list_of_files)
    i = load_untouch_nii(list_of_files(j).name);
    I = i.img;
    I_M = I.*mask;
    vt = reshape(I_M, [1,xres*yres*zres]);
    [~,idx] = find(vt);
    vt = vt(:,idx);
    data(j,:) = vt;
end
% Plotting the ks density curves to identify the best subset of patients
cd('F:\LRCBH\Projects\COBRE\Results\Documents\Excel');
%fluency_ratio = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 3 107 3]);
% pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 9 90 9]);
% pf_avg = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 10 90 10]);
% sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v2','Range',[2 11 90 11]);
pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v3','Range',[2 9 44 9]);
pf_avg = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v3','Range',[2 10 44 10]);
sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression_v3','Range',[2 11 44 11]);
[pf1, xi1] = ksdensity(pf(1:26, :));
[pf2, xi2] = ksdensity(pf(27:44, :));
plot(xi1,pf1,xi2,pf2)
lgd = legend('Control','MCI');
title(lgd, 'Phonemic fluency')
figure;
[sf1, xi3] = ksdensity(sf(1:26, :));
[sf2, xi4] = ksdensity(sf(27:43, :));
plot(xi3,sf1,xi4,sf2)
lgd = legend('Control','MCI');
title(lgd, 'Semantic fluency')
figure;
[pf_avg1, xi5] = ksdensity(pf_avg(1:26, :));
[pf_avg2, xi6] = ksdensity(pf_avg(27:43, :));
plot(xi5,pf_avg1,xi6,pf_avg2)
lgd = legend('Control','MCI');
title(lgd, 'Mean Phonemic fluency')
% Group specific maps
figure;
[sf_grp1, xi7] = ksdensity(sf(1:44, :));
[pf_grp1, xi8] = ksdensity(pf_avg(1:44, :));
plot(xi7,sf_grp1,xi8,pf_grp1)
lgd = legend('Semantic','Phonemic');
title(lgd, 'Controls')
figure;
[sf_grp2, xi9] = ksdensity(sf(45:89, :));
[pf_grp2, xi10] = ksdensity(pf_avg(45:89, :));
plot(xi9,sf_grp2,xi10,pf_grp2)
lgd = legend('Semantic','Phonemic');
title(lgd, 'MCI')

% Testing the regression model
mdl = fitlm(X, Y(:, 2));
