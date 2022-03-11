clc
clear
tic

%------------------------------------------------------------------------------%
%                     IMPORTING LIST OF USEFUL DIRECTORIES
%------------------------------------------------------------------------------%

rootdir = 'E:\LRCBH\Data\COBRE-MNI\Individual_data\Useful';
pca_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\1.PCA\Useful';
ica_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\1.ICA';
dr_savedir = 'E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\2.DR';
fcn_savedir='E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\3.FCN';
asso_savedir = 'E:\LRCBH\Results\Matlab\v2\5.Association';
dirpath = dir(rootdir);
subdir = [dirpath(:).isdir];
subloc = {dirpath(subdir).name}';
subloc(ismember(subloc,{'.','..'})) = [];

%------------------------------------------------------------------------------%
%                     ACQUIRING CRUCIAL VOXEL INFORMATION
%------------------------------------------------------------------------------%

cd('E:\LRCBH\Data\MNI_segmented');
%cd('W:\MNI_segmented')
m = load_untouch_nii('standard_binary.nii');
M = m.img;
nii_temp = m;
[x, y, z] = size(M);
mask_temp = reshape(M, [1,x*y*z]);
[~, indices] = find(mask_temp);
trs = 850;
comp = 200;

%------------------------------------------------------------------------------%
%                                 PCA REDUCTION
%------------------------------------------------------------------------------%

%Iterating through each of the subjects
for i=1:length(subloc)
    suboi = subloc{i};
    current = strcat(rootdir,'\', suboi);
    cd(current);
    files = dir('*.nii');
    temp = zeros(trs, 228453);
    fprintf('Generating voxeltime data for subject %s...\n',suboi);
    [temp,index] = generate_vt(files,M,x,y,z);
    vt_data = double(zscore(temp(16:end, :)));
    clear temp;
    fprintf('Performing PCA reduction using %d components for subject %s...\n',comp,suboi);
    PCA_red = double(do_PCA(vt_data,comp));
    cd(pca_savedir);
    save(fullfile(pca_savedir, sprintf('PCA_%s.mat',suboi)),'PCA_red','vt_data','index');
end
fprintf('PCA reduction done.\n');

%------------------------------------------------------------------------------%
%                               TEMPORAL CONCATENATION
%------------------------------------------------------------------------------%

fprintf('Performing temporal concatenation of subjects...\n');
tcat_data = temporal_concat(subloc,pca_savedir);
fprintf('Temporal concatenation done.\n')
clearvars -except tcat_data ica_savedir

%------------------------------------------------------------------------------%
%                               SPATIAL ICA BASED ON HTANH
%------------------------------------------------------------------------------%

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
save(fullfile(ica_savedir,sprintf('gica_%d_result.mat',npca)),'A','S','W','E','eigval','White','-v7.3');
toc

%------------------------------------------------------------------------------%
%                  GENERATING NII VERSION OF ICA COMPONENTS
%------------------------------------------------------------------------------%

cd(ica_savedir);
ica_dat = load('gICA_100_result.mat','S');
S = ica_dat.S;
S = double(S);
save_ica_nii(S,x,y,z,indices,m,'gICA',ica_savedir);

%------------------------------------------------------------------------------%
%                                      DUAL REGRESSION
%------------------------------------------------------------------------------%

%noise_idx = [2,4,11,13,14,17,18,19,28,23];
%Using only the useful components
%S(noise_idx,:)=[];
%Getting the directory having the subject wise voxel time data
dirloc = dir(pca_savedir);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for i =1:length(subloc)
    suboi = subloc{i};
    current = strcat(pca_savedir,'\', suboi);
    sub_data = load(current,'vt_data');
    sub_data = (sub_data.vt_data)';
    fprintf('Performing dual regression for subject %s...\n',suboi(5:7));
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
    save_ica_nii(ss_sm',x,y,z,indices,m,'ica_',ss_dr_dir);
end
fprintf('Dual regression done.\n')

%------------------------------------------------------------------------------%
%                          FUNCTIONAL CONNECTIVITY ANALYSIS
%------------------------------------------------------------------------------%

dirpath = dir(dr_savedir);
subdir = [dirpath(:).isdir];
subpath = {dirpath(subdir).name}';
subpath(ismember(subpath,{'.','..'})) = [];
for i=1:length(subpath)
    sub = subpath{i};
    current = strcat(dr_savedir,'\', sub);
    cd(current);
    dr = load('dualregression.mat','ss_tc');
    tc = (dr.ss_tc)';
    [pval, fcn] = corrcoef(tc);
    save(fullfile(fcn_savedir, sprintf('FCN_%s.mat',sub)),'fcn','pval');
end

%------------------------------------------------------------------------------%
%                          CORRELATION WITH FLUENCY SCORE
%------------------------------------------------------------------------------%

fluency_dir = 'E:\LRCBH\Projects\COBRE\Results\Documents\Excel';
% Specifying the corresponding patient indices
noise_idx =[2,12,17,25,28,31,32,35,44,47,51,52,55,62,66,67,71,78,79,84,85,93,96,99]; % ICA 100
noise_idx = [];
all_idx = [1:102]';
grp1_idx = [1,2,5,9,12,17,20,22,30,31,32,36,38,39,42,43,44,45,46,48,49,50,51,54,56,57,58,...
    59,60,61,62,64,67,69,72,75,76,78,80,81,82,83,84,86,87,88,89,94,100,101,102];
grp2_idx = setdiff(all_idx,grp1_idx)';
% Declaring other relevant information
edges = -0.6 : 0.2 : 0.6;
pos_range = [0.4, 0.6];
neg_range = [-0.6, -0.4];
% Getting the results
[grp1.fluency_ratio, grp1.subinfo, grp1.fcn_matrix, grp1.rsn_loc, grp1.correlation_coefficient, ...
    grp1.min_correlation, grp1.max_correlation, grp1.number, grp1.poscorrelation_value, grp1.negcorrelation_value, grp1.pos_rsn, grp1.neg_rsn] = ...
    fcnmat_results(fcn_savedir, fluency_dir, grp1_idx,noise_idx,1, edges, pos_range, neg_range);
[grp2.fluency_ratio, grp2.subinfo, grp2.fcn_matrix, grp2.rsn_loc, grp2.correlation_coefficient, ...
    grp2.min_correlation, grp2.max_correlation, grp2.number, grp2.poscorrelation_value, grp2.negcorrelation_value, grp2.pos_rsn, grp2.neg_rsn] = ...
    fcnmat_results(fcn_savedir, fluency_dir, grp2_idx,noise_idx,2, edges, pos_range, neg_range);
[grp_all.fluency_ratio, grp_all.subinfo, grp_all.fcn_matrix, grp_all.rsn_loc, grp_all.correlation_coefficient, ...
    grp_all.min_correlation, grp_all.max_correlation, grp_all.number, grp_all.poscorrelation_value, grp_all.negcorrelation_value, grp_all.pos_rsn, grp_all.neg_rsn] = ...
    fcnmat_results(fcn_savedir, fluency_dir, all_idx,noise_idx,3, edges, pos_range, neg_range);
save(fullfile('E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results',sprintf('Correlation_fluency_results_%d.mat',100)),'grp1','grp2','grp_all');

%------------------------------------------------------------------------------%
%   BAR PLOT SHOWING THE DISTRIBUTION OF CORRELATION VALUES
%------------------------------------------------------------------------------%

xloc = [1,2,3,4,5,6];
% Grp1
grp1_number = grp1.number;
grp1_percent = (grp1_number/sum(grp1_number))*100;
figure;
bar(xloc, grp1_percent, 0.5);
ax=gca;
ax.FontSize=24;
xticklabels({'[-0.6 : -0.4]','[-0.4 : -0.2]','[-0.2 : 0.0]','[0.0 : 0.2]','[0.2 : 0.4]','[0.4 : 0.6]'});
ylabel(' Percentage of connections ','fontweight','bold','FontSize',32);
xlabel(' Correlation coefficient : Range ','fontweight','bold','FontSize',32);
title(' Grp I - Normal Cognition','fontweight','bold','FontSize',40);
% Grp2
grp2_number = grp2.number;
grp2_percent = (grp2_number/sum(grp2_number))*100;
figure;
bar(xloc, grp2_percent, 0.5);
ax=gca;
ax.FontSize=24;
xticklabels({'[-0.6 : -0.4]','[-0.4 : -0.2]','[-0.2 : 0.0]','[0.0 : 0.2]','[0.2 : 0.4]','[0.4 : 0.6]'});
ylabel(' Percentage of connections ','fontweight','bold','FontSize',32);
xlabel(' Correlation coefficient : Range ','fontweight','bold','FontSize',32);
title(' Grp II - MCI','fontweight','bold','FontSize',40);

%---------------------Boxplot with data scattared-------------------------%
a1 = (grp1.correlation_coefficient)';
a2 = (grp2.correlation_coefficient)';
a = vertcat(a1,a2);
g1 = [zeros(length(a1),1); ones(length(a2),1)];
boxplot(a,g1,'Notch','off');
%set(gca,'YLim',[-3 3])
set(gca,'XTickLabel', {'Normal cognition','MCI'});
set(gca,'TickLabelInterpreter','none');
set(gca,'fontsize',20);
hold on
[C,~,ic] = unique(g1,'stable');
scatter(ic,a,'k','filled','MarkerFaceAlpha',0.5');
color = ['g','r'];
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.2);
 end
ylabel('Correlation coefficient','fontweight','bold','FontSize',24)
title('f2','fontweight','bold','FontSize',24)
legend(h([3,4]),{'Diseased','Healthy'})
% Generating the barplots with bins denoting eh number of connections having
% that correlation value


%------------------------------------------------------------------------------%
%                                REGRESSION ANALYSIS
%------------------------------------------------------------------------------%

% Generating the design matrix
fluency_ratio = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 3 103 3]);
grp = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 4 103 4]);
age = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 6 103 6]);
ed = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 7 103 7]);
suvr = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 8 103 8]);
pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 9 103 9]);
sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 10 103 10]);
suvr_dis = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 11 103 11]);

%regressor = horzcat(fluency_ratio, grp);
regressor = horzcat(pf, sf, grp);
covariates = horzcat(age, ed, suvr_dis);
interaction = fluency_ratio.*grp;

% Fitting regression model for all the components
rsn_idx = [17:20];
for i=1:length(rsn_idx)
    idx = rsn_idx(i);
    fprintf('Fitting regression model for component %d...\n',idx);
    [~,~,~,~,~,~,~,~,~,~,~,~,~,~] = regress_model('E:\LRCBH\Results\Matlab\3.DR\Unbiased', ...
        regressor, [], covariates,0.05,'dualregression',idx,'E:\LRCBH\Results\Matlab\v2\8.Association_separate\Regression_models');
end
fprintf('Fitting regression model done...\n');

%x = 1-normcdf(5);

% Testing the assumptions of the regression model
% 1. Linear relationship
Y_mean = mean(Y,2);
X = horzcat(regressor, interaction, covariates);
for j = 1:size(X,2)
    figure;
    scatter(X(:,j),Y_mean);
end
scatter(X(:,6),Y_mean,'filled');
lsline;
title('Testing linear relationship');
ylabel('y - Dependent variable')
%xlabel('x1 - Fluency ratio');
%xlabel('x2 - Group ID');
%xlabel('x3 - Interaction term');
%xlabel('x4 - Age');
%xlabel('x5 - Education');
xlabel('x6 - SUVR');
x2 = horzcat(Y_mean(1:44,:),Y_mean(45:end,:));
boxplot(x2);
x6 = sqrt(X(:,6));
scatter(x6,Y_mean,'filled');

% Generating scatter plots across groups for significant voxels
Y_sig = Yfitted(:,sig_voxel(1,:));
Y_avg = mean(Y_sig,2);
gscatter(fluency_ratio,Y_avg,grp);
lsline;
xlabel('Fluency ratio');
ylabel('Mean voxel value')
legend('MCI','Normal');
% Generating mean spatial maps across groups
Y_sig_grp = mean(Y_sig);
loc = sig_voxel(1,:);
temp = zeros(1,228453);
temp(1,loc) = Y_sig_grp;
%Saving them as nifti files
save_ica_nii(temp,x,y,z,indices,m,'association_grp','E:\LRCBH\Results\Matlab\v2\5.Association\Spatial_maps');
% Plotting the scatter plot
Y_sig = regress_withinter.Yfitted(:,sig_voxels(1,:));
Y_mean = mean(Y_sig,2);
gscatter(fluency_ratio,Y_mean,grp);
lsline;
xlabel('Fluency ratio');
ylabel('Mean voxel value')
legend('MCI','Normal');
pval_interaction = (regress_withinter.pvalue(:,2))';% For the interaction term
%alpha_level = 0.05/n_compar;
sig_idx= find(pval_interaction<0.05);
%Testing the linear relationship of the residuals vs fitted
scatter(mean(regress_withinter.Yfitted,2),mean(regress_withinter.residuals,2));
lsline;
%Q-Q plot of the residuals
qqplot(mean(regress_withinter.residuals,2));
msqe = mse(regress_withinter.Yfitted,Y);
scatter(mean(regress_withinter.Yfitted,2))

% Finding the significant voxels based on alpha level and the Rsquared and
% Fsquared statistic of those components
loc = 'E:\LRCBH\Results\Matlab\v2\5.Association\Regression_models';
dirloc = dir(loc);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
pval_voxel=  0.001;
saveloc = 'E:\LRCBH\Results\Matlab\v2\5.Association\Spatial_maps';
for j = 1:length(subloc)
    suboi = subloc{j};
    current = strcat(loc,'\',suboi);
    regress_result = load(current);
    % Finding the significant voxels with grp
    idx_grp = find(regress_result.pval(:,3)<0.001)';
    vox_grp{j,1} = idx_grp;
    % Generating t-statistic map - group significance
    [~,~,~,stats] = ttest2(regress_result.Yfitted(1:44,idx_grp),regress_result.Yfitted(45:end,idx_grp));
    temp = zeros(1,228453);
    temp(1,idx_grp) = stats.tstat;
    save_ica_nii(temp,x,y,z,indices,m,'tstat_map',strcat(saveloc,'\','Group','\',suboi(21:22)));
    clear temp stats;
    % Finding the significant voxels for interaction
    idx_interaction = find(regress_result.pval(:,4)<0.001)';
    vox_interaction{j,1} = idx_interaction;
    % Generating t-statistic map - interaction significance
    [~,~,~,stats] = ttest2(regress_result.Yfitted(1:44,idx_interaction),regress_result.Yfitted(45:end,idx_interaction));
    temp = zeros(1,228453);
    temp(1,idx_interaction) = stats.tstat;
    save_ica_nii(temp,x,y,z,indices,m,'tstat_map',strcat(saveloc,'\','Interaction','\',suboi(21:22)));
    clear temp stats;
    % Finding the significant voxels common across grp & interaction
    idx_common = intersect(idx_grp,idx_interaction)';
    vox_common{j,1} = idx_common;
    % Generating t-statistic map - common significance
    [~,~,~,stats] = ttest2(regress_result.Yfitted(1:44,idx_common),regress_result.Yfitted(45:end,idx_common));
    temp = zeros(1,228453);
    temp(1,idx_common) = stats.tstat;
    save_ica_nii(temp,x,y,z,indices,m,'tstat_map',strcat(saveloc,'\','Common','\',suboi(21:22)));
    clear temp stats;
    % Finding the effect size
    R_grp = (regress_result.Rsquared_adjust(idx_grp,1))';
    R_interaction = (regress_result.Rsquared_adjust(idx_interaction,1))';
    R_common = (regress_result.Rsquared_adjust(idx_common,1))';
    R_grp_mean = mean(R_grp,2);
    R_interaction_mean= mean(R_interaction,2);
    R_common_mean= mean(R_common,2);
    f_grp(j,1)= R_grp_mean/(1-R_grp_mean);
    f_interaction(j,1) = R_interaction_mean/(1-R_interaction_mean);
    f_common(j,1) = R_common_mean/(1-R_common_mean);
end
save(fullfile('E:\LRCBH\Results\Matlab\v2\5.Association', sprintf('mlr_model_results_all.mat')),'subloc','vox_common','vox_grp','vox_interaction','f_grp','f_interaction','f_common');
for i = 1:size(vox_grp,1)
    vox_grp_n(i,1) = numel(vox_grp{i,1});
    vox_interaction_n(i,1) = numel(vox_interaction{i,1});
    vox_common_n(i,1) = numel(vox_common{i,1});
end

% Generating scatter plots for the ICA components
loc = 'E:\LRCBH\Results\Matlab\v2\5.Association\Regression_models';
dirloc = dir(loc);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
subloi = subloc{3,1};
cd('E:\LRCBH\Results\Matlab\v2\5.Association\Scatter_plots\Common');
for j = 1:size(subloc,1)
    suboi = subloc{j};
    current = strcat(loc,'\',suboi);
    regress_result = load(current);
    Y_sig = regress_result.Yfitted(:,vox_common{j,1});
    Y_mean = mean(Y_sig,2);
    figure;
    gscatter(fluency_ratio,Y_mean,grp);
    lsline;
    xlabel('Fluency ratio');
    ylabel('Mean voxel value')
    legend('Normal','MCI');
    title(sprintf('Component %s',suboi(21:22)));
    %saveas(h,sprintf('Component %s'),suboi(21:22));
    clear Y_sig Y_mean suboi;
end
% Generate ICA maps based on t-statistics of the interaction term
loc = 'E:\LRCBH\Results\Matlab\v2\6.Association_v2\Regression_models';
dirloc = dir(loc);
subloc = {dirloc.name}';
subloc(ismember(subloc,{'.','..'})) = [];
for k = 1:length(subloc)
    suboi = subloc{k};
    current = strcat(loc,'\',suboi);
    regress_result = load(current);
    Tstat(k,:) = regress_result.tstat(:,1)';
end
save_ica_nii(Tstat,x,y,z,indices,m,'tmap_interaction_','E:\LRCBH\Results\Matlab\v2\6.Association_v2\Spatial maps');



% Differences in the distribution of the fluency score
fluency_ratio = readmatrix('Cobre_fluency_study.xlsx','Sheet','regression','Range',[2 3 89 3]);
% pf = readmatrix('Cobre_language_study.xlsx','Sheet','regression','Range',[2 8 89 8]);
pf = readmatrix('Cobre_fluency_study.xlsx','Sheet','regression','Range',[2 9 89 9]);
sf = readmatrix('Cobre_fluency_study.xlsx','Sheet','regression','Range',[2 10 89 10]);
[h,p,ci,stats] = ttest2(fluency_ratio(1:44,1),fluency_ratio(45:end,1));
[h1,p1,ci1,stats1] = ttest2(pf(1:44,1),pf(45:end,1));
[h2,p2,ci2,stats2] = ttest2(sf(1:44,1),sf(45:end,1));


%ks density plots
% Phonemic Fluency
[f,xi] = ksdensity(pf);
figure
plot(xi,f);
title('Phonemic fluency');
clear f xi;
% Semantic Fluency
[f,xi] = ksdensity(sf);
figure
plot(xi,f);
title('Semantic fluency');
clear f xi;
% Fluency ratio
[f,xi] = ksdensity(fluency_ratio);
figure
plot(xi,f);
title('Fluency ratio');
clear f xi;
% Group
[f,xi] = ksdensity(grp);
figure
plot(xi,f);
title('Group');
clear f xi;
% Age
[f,xi] = ksdensity(age);
figure
plot(xi,f);
title('Age');
clear f xi;
% Education
[f,xi] = ksdensity(ed);
figure
plot(xi,f);
title('Education');
clear f xi;
% SUVR - continuous
[f,xi] = ksdensity(suvr);
figure
plot(xi,f);
title('SUVR - continuous');
clear f xi;
% SUVR - discrete
[f,xi] = ksdensity(suvr_dis);
figure
plot(xi,f);
title('SUVR - discrete');
clear f xi;

% Checking the correlation coefficient and VIF
mat = horzcat(fluency_ratio, interaction, age, ed, suvr_dis);
mat = horzcat(regressor, covariates);
R= corrcoef(mat);
V=diag(inv(R))';

% Calculating the p-value
[h,p] = ttest(S(83,:),S(50,:));


% Testing out the regression analysis
