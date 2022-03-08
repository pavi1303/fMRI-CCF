function oobErr = oobErrRF(params,X)
randomForest = TreeBagger(100,X,'MPG','Method','regression',...
    'OOBPrediction','on','MinLeafSize',params.minLS);
oobErr = oobQuantileError(randomForest);
end