%------------------------------------------------------------------------------%
%  IDENTIFICATION OF THE APPROPRIATE TRANSFORMATION METHOD
%------------------------------------------------------------------------------%
%------------------------------------------------------------------------------%
%                                    SEMANTIC FLUENCY
%------------------------------------------------------------------------------%
% No transformation
cd('E:\LRCBH\Projects\COBRE\Results\Documents\Excel');
sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 10 107 10]);
% Log transformation
sf_l = log(sf);
% Box-cox transformation
[sf_bc, lambda_sf] = boxcox(sf);
% Semantic fluency - Generating the ks density functions
% No transform
[f1, xi1] = ksdensity(sf(1:51,1));
[f2, xi2] = ksdensity(sf(52:end,1));
figure;
subplot(2,2,1);
plot(xi1,f1,xi2,f2);
lgd = legend('NC','MCI');
title('Semantic fluency (No transform)');
% Log transform
[f1, xi1] = ksdensity(sf_l(1:51,1));
[f2, xi2] = ksdensity(sf_l(52:end,1));
subplot(2,2,2);
plot(xi1,f1,xi2,f2);
lgd = legend('NC','MCI');
title('Semantic fluency (Log transform)');
% Box-cox transform
[f1, xi1] = ksdensity(sf_bc(1:51,1));
[f2, xi2] = ksdensity(sf_bc(52:end,1));
subplot(2,2,3);
plot(xi1,f1,xi2,f2);
lgd = legend('NC','MCI');
title('Semantic fluency (Box-cox transform)');

% Natural log transformation
Y_l = log(Y);
% Scatter plot
%1. No transformation
figure;
gscatter(sfl,Y(:,1),X(:,2));
title('Level-log transform');
ylabel('Mean FC value (natural log)');
xlabel('Semantic fluency');
% Getting the untransformed semantic fluency
cd('E:\LRCBH\Projects\COBRE\Results\Documents\Excel');
sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 10 107 10]);
sfl = log(sf);
% Level-log transformation
figure;
subplot(2,2,1);
gscatter(sfl,Y(:,1),X(:,2));
legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
title('Level-log transform');
ylabel('Mean FC value');
xlabel('Semantic fluency (Natural log)');
% Log-level transform
subplot(2,2,2);
gscatter(sf,Y_l(:,1),X(:,2));
legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
title('Log-level transform');
ylabel('Mean FC value (Natural log)');
xlabel('Semantic fluency ');
% Log-log transform
subplot(2,2,3);
gscatter(sfl,Y_l(:,1),X(:,2));
legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
title('Log-log transform');
ylabel('Mean FC value (Natural log)');
xlabel('Semantic fluency (Natural log)');

%------------------------------------------------------------------------------%
%                                    TARGET VARIABLE
%------------------------------------------------------------------------------%
% Log transform
Y_l = log(Y);
% Box-cox transform
for i=1:size(Y,2)
    [Y_bc(:,i), lambda(:,i)] = boxcox(Y(:,i));
end
cd('E:\LRCBH\Weekly_meeting\Images\0329\Target variable');
for rsn_idx = 1:14
    h=figure;
    subplot(2,2,1);
    [f3, xi3] = ksdensity(Y(1:51,rsn_idx));
    [f4, xi4] = ksdensity(Y(52:end,rsn_idx));
    plot(xi3,f3,xi4,f4);
    title('No transformation');
    subplot(2,2,2);
    [f3, xi3] = ksdensity(Y_l(1:51,rsn_idx));
    [f4, xi4] = ksdensity(Y_l(52:end,rsn_idx));
    plot(xi3,f3,xi4,f4);
    title('Log transformation');
    subplot(2,2,3);
    [f3, xi3] = ksdensity(Y_bc(1:51,rsn_idx));
    [f4, xi4] = ksdensity(Y_bc(52:end,rsn_idx));
    plot(xi3,f3,xi4,f4);
    title('Box cox transformation');
    %saveas(h,sprintf('RSN_%d.jpg',rsn_idx));
end
%------------------------------------------------------------------------------%
%                                    SCATTER PLOTS
%------------------------------------------------------------------------------%
cd('E:\LRCBH\Weekly_meeting\Images\0329\Scatter plots');
for rsn_idx = 1:14
    figure;
    subplot(2,2,1);
    gscatter(X(:,1),Y_l(:,rsn_idx),X(:,2));
    legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
    title('Log-level transform');
    ylabel('Mean FC (Natural log)');
    xlabel('Semantic fluency ');
    % Log-log transform
    subplot(2,2,2);
    gscatter(X_l(:,1),Y_l(:,rsn_idx),X(:,2));
    legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
    title('Log-log transform');
    ylabel('Mean FC (Natural log)');
    xlabel('Semantic fluency (Natural log)');
    % Boxcox-level transform
    subplot(2,2,3);
    gscatter(X(:,1),Y_bc(:,rsn_idx),X(:,2));
    legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
    title('Boxcox-level transform');
    ylabel('Mean FC (Box-cox)');
    xlabel('Semantic fluency ');
    % Boxcox-boxcox transform
    subplot(2,2,4);
    gscatter(X_l(:,1),Y_bc(:,rsn_idx),X(:,2));
    legend({'Normal Cognition', 'MCI'}, 'Location','southoutside');
    title('Boxcox-boxcox transform');
    ylabel('Mean FC (Box-cox)');
    xlabel('Semantic fluency (Box-cox)');
    %saveas(h,sprintf('Scatterplot_%d.jpg',rsn_idx));
end
% Original X variable
cd('E:\LRCBH\Projects\COBRE\Results\Documents\Excel');
fluency_ratio = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 3 107 3]);
grp = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 4 107 4]);
age = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 6 107 6]);
ed = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 7 107 7]);
%interaction = readmatrix('Cobre_fluency_study_v2.xlsx', 'Sheet','regression','Range',[2 5 107 5]);
%suvr = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 8 103 8]);
pf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 9 107 9]);
sf = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 10 107 10]);
suvr_dis = readmatrix('Cobre_fluency_study_v2.xlsx','Sheet','regression','Range',[2 11 107 11]);

%regressor = horzcat(fluency_ratio, grp);
%regressor = horzcat(fluency_ratio, grp, interaction);
interaction = grp.*sf;
regressor = horzcat(sf, grp, interaction);
covariates = horzcat(pf,age, ed, suvr_dis);
X = horzcat(regressor, covariates);
% X - Log transformation
X_l = zeros(106,7);
X_l(:,1) = log(X(:,1));
X_l(:,2) = X(:,2);
X_l(:,3) = X_l(:,1).*X_l(:,2);
X_l(:,4) = log(X(:,4));
X_l(:,5) = log(X(:,5));
X_l(:,6) = log(X(:,6));
X_l(:,7) = X(:,7);
% X - Box cox transformation
X_bc = zeros(106,7);
[X_bc(:,1), lambda] = boxcox(X(:,1));
X_bc(:,2) = X(:,2);
X_bc(:,3) = X_bc(:,1).*X_bc(:,2);
[X_bc(:,4), lambda] = boxcox(X(:,4));
[X_bc(:,5), lambda]= boxcox(X(:,5));
[X_bc(:,6), lambda] = boxcox(X(:,6));
X_bc(:,7) = X(:,7);
inputdata_savedir = cd('E:\LRCBH\Projects\COBRE\Results\Matlab\ICA_100_results\5.Regression\ROI_analysis');
save(fullfile(inputdata_savedir, sprintf('ROI_regression_input_data.mat')),'X','Y','X_l','Y_l','X_bc','Y_bc');

Q = quantile(X(:,1),[0.25 0.5 0.75]);
% Trying cubic smoothing spline
spline = csaps(X(1:51,1),Y_l(1:51,1),1);
spline1 = csaps(X(52:end,1),Y_l(52:end,1),1);
figure;
fnplt(spline);
hold on
fnplt(spline1);
plot(X(1:51,1),Y_l(1:51,1),'ko');
plot(X(52:end,1),Y_l(52:end,1),'ro');
hold off
title('Cubic smoothing spline - Trial');
X1_nc = sort(X(1:51,1));
X1_mci = sort(X(52:end,1));

