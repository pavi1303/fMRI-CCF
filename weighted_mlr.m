for j = 1:length(subloc)
    suboi = subloc{j};
    current = strcat(loc,'\',suboi);
    regress_result = load(current);
    % Finding the significant voxels with grp
    idx_grp = find(regress_result.pval(:,3)<0.001)';
    %vox_grp{j,1} = idx_grp;
    % Generating t-statistic map - group significance
    [~,~,~,stats] = ttest2(regress_result.Yfitted(1:44,idx_grp),regress_result.Yfitted(45:end,idx_grp));
    temp = zeros(1,228453);
    temp(1,idx_grp) = stats.tstat;
    save_ica_nii(temp,x,y,z,indices,m,'tstat_map',strcat(saveloc,'\','Group','\',suboi(21:22)));
    clear temp stats;
    % Finding the significant voxels for interaction
    idx_interaction = find(regress_result.pval(:,4)<0.001)';
    %vox_interaction{j,1} = idx_interaction;
    % Generating t-statistic map - interaction significance
    [~,~,~,stats] = ttest2(regress_result.Yfitted(1:44,idx_interaction),regress_result.Yfitted(45:end,idx_interaction));
    temp = zeros(1,228453);
    temp(1,idx_interaction) = stats.tstat;
    save_ica_nii(temp,x,y,z,indices,m,'tstat_map',strcat(saveloc,'\','Interaction','\',suboi(21:22)));
    clear temp stats;
    % Finding the significant voxels common across grp & interaction
    idx_common = intersect(idx_grp,idx_interaction)';
    %vox_common{j,1} = idx_common;
    % Generating t-statistic map - common significance
    [~,~,~,stats] = ttest2(regress_result.Yfitted(1:44,idx_common),regress_result.Yfitted(45:end,idx_common));
    temp = zeros(1,228453);
    temp(1,idx_common) = stats.tstat;
    save_ica_nii(temp,x,y,z,indices,m,'tstat_map',strcat(saveloc,'\','Interaction','\',suboi(21:22)));
    clear temp stats;
    % Finding the effect size
%     R_grp = (regress_result.Rsquared_adjust(idx_grp,1))';
%     R_interaction = (regress_result.Rsquared_adjust(idx_interaction,1))';
%     R_common = (regress_result.Rsquared_adjust(idx_common,1))';
%     R_grp_mean = mean(R_grp,2);
%     R_interaction_mean= mean(R_interaction,2);
%     R_common_mean= mean(R_common,2);
%     f_grp(j,1)= R_grp_mean/(1-R_grp_mean);
%     f_interaction(j,1) = R_interaction_mean/(1-R_interaction_mean);
%     f_common(j,1) = R_common_mean/(1-R_common_mean);
end