function [acc, dcc, u, l] = accordance_discordance_estimate(ts_data, q)
% STEP 0 : Checking the accordance discordance estimation

% STEP 1: Time series normalization
for ts = 1:size(ts_data, 1)
    ts_norm(ts, :) = zscore(ts_data(ts, :));
end
% STEP 2: Identification of the extreme points in the time series
u = norminv(q);
l = norminv(1-q);
% STEP 3: Activation vector formation
for i = 1:size(ts_norm, 1)
    trs = ts_norm(i, :)
    for j = 1: length(trs)
        if trs(1, j)<u
            xu(i, j) = 0;
        else
            xu(i, j) = 1;
        end
    end
end
% STEP 4: Deactivation vector formation
for i = 1:size(ts_norm, 1)
    trs = ts_norm(i, :)
    for j = 1: length(trs)
        if trs(1, j)>l
            xl(i, j) = 0;
        else
            xl(i, j) = -1;
        end
    end
end
% STEP 5: Calculation of accordance and discordance values
for i = 1:size(ts_norm, 1)
    for j = 1:size(ts_norm, 1)
        E(i, j) = (sqrt((xu(i, :)*xu(i, :)' + xl(i, :)*xl(i, :)'))*sqrt((xu(j, :)*xu(j, :)' + xl(j, :)*xl(j, :)')));
        acc(i, j) = (xu(i, :)*xu(j, :)' + xl(i, :)*xl(j, :)')/E(i, j);
        dcc(i, j) = (xu(i, :)*xl(j, :)' + xl(i, :)*xu(j, :)')/E(i, j);
    end
end
for i = 1:length(acc)
    acc(i, i) = round(acc(i, i));
end
% STEP 6 : Checking the estimated accordance and discordance matrices
if unique(diag(acc)) == 1 & unique(diag(dcc)) == 0
    fprintf('Accordance and discordance estimation ran succesfully!!!\n');
else
    warning('Error in estimation of the accordance and discordance values');
end
end
