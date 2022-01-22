function [data_concat] = temporal_concat(files_list, savepath)
for i=1:length(files_list)
    pat = files_list{i};
    file = fullfile(savepath, sprintf('PCA_%s.mat', pat));
    data = load(file);
    data = struct2cell(data);
    pat_data{1,i} = data;
end
pca_tcat = cellfun(@vertcat, pat_data);
data_concat = vertcat(pca_tcat{:});
end
