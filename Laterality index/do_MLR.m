function mlr_results = do_MLR(cov_mat, grp, sf, pf, fc_mat, fluency_var, n_ica, noise_idx, varoi_l1, varoi_l2)
mlr_results = struct;
if fluency_var=='SF'
    Y =sf;
elseif fluency_var=='PF'
    Y =pf;
end
% Type 1 - Traditional regression model Var
mlr_level1 = struct;
for j = 1:size(fc_mat, 2)
    X = [cov_mat, grp, fc_mat(:, j), grp.*fc_mat(:,j)];
    mlr_level1.X = X;
    mlr_level1.Y = Y;
    %mdl = fitlm(X, Y, 'ResponseVar', {'SF'}, 'PredictorVars','PredictorVars', {'Age', 'Education', 'Gender', 'Group', 'FC', 'Interaction'}, 'CategoricalVar',{'Group'});
    mdl = fitlm(X, Y);
    mlr_level1.adj_r2(j,1) = mdl.Rsquared.Adjusted;
    pvalue(j,:) = mdl.Coefficients.pValue;
    tstat(j,:) = mdl.Coefficients.tStat;
    mlr_level1.beta(j,:) = mdl.Coefficients.Estimate;
    mlr_level1.Yfitted(:,j) = mdl.Fitted;
    mlr_level1.raw_res(:,j) = mdl.Residuals.Raw;
end
mlr_level1.pvalue = array2table(pvalue,'VariableNames', [{'Intercept', 'Age', 'Education', 'Gender', 'Group', 'FC', 'Interaction'}]);
mlr_level1.Tstat = array2table(tstat,'VariableNames', [{'Intercept', 'Age', 'Education', 'Gender', 'Group', 'FC', 'Interaction'}]);
% Getting the indices of the significant p-values if any
new_ica = setdiff([1:30], noise_idx);
[loc(:,1), loc(:, 2)] = find(tril(rand(length(new_ica)), -1));
mlr_level1.sig_connection = new_ica(loc(find(pvalue(:, varoi_l1)<0.05), :));
mlr_level1.sig_pvalue = array2table(pvalue(find(pvalue(:, varoi_l1)<0.05), :),...
    'VariableNames', [{'Intercept', 'Age', 'Education', 'Gender', 'Group', 'FC', 'Interaction'}]);

clearvars pvalue tstat mdl;
% Type 2 - Two level regression model
cov_mdl = fitlm(cov_mat, Y); % Regressing out the effects of covariates
Y1 = cov_mdl.Residuals.Raw;
mlr_level2 = struct;
for j = 1:size(fc_mat, 2)
    X = [grp, fc_mat(:, j), grp.*fc_mat(:,j)];
    mlr_level2.X = X;
    mlr_level1.Y = Y1;
    mdl = fitlm(X, Y1, 'ResponseVar',fluency_var, 'PredictorVars',{'Group', 'FC', 'Interaction'});
    mlr_level2.adj_r2(j,1) = mdl.Rsquared.Adjusted;
    pvalue(j,:) = mdl.Coefficients.pValue;
    tstat(j,:) = mdl.Coefficients.tStat;
    mlr_level2.beta(j,:) = mdl.Coefficients.Estimate;
    mlr_level2.Yfitted(:,j) = mdl.Fitted;
    mlr_level2.raw_res(:,j) = mdl.Residuals.Raw;
end
mlr_level2.pvalue = array2table(pvalue,'VariableNames', [{'Intercept'}, {'Group', 'FC', 'Interaction'}]);
mlr_level2.Tstat = array2table(tstat,'VariableNames', [{'Intercept'}, {'Group', 'FC', 'Interaction'}]);
% Getting the indices of the significant p-values if any
mlr_level2.sig_connection = new_ica(loc(find(pvalue(:, varoi_l2)<0.05), :));
mlr_level2.sig_pvalue = array2table(pvalue(find(pvalue(:, varoi_l2)<0.05), :),...
    'VariableNames', [{'Intercept', 'Group', 'FC', 'Interaction'}]);

% The output structure with results for both types of models
mlr_results.type1 = mlr_level1;
mlr_results.type2 = mlr_level2;
end