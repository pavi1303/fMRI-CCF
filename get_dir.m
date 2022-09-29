function [folders] = get_dir(location)
cd(location);
contents = dir(location);
temp = [contents(:).isdir];
folders = {contents(temp).name}';
folders(ismember(folders,{'.','..'})) = [];
end


