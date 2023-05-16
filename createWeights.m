function Weights = createWeights(var,vars,addVars,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

tic

climShort = strcat(vars,'_ShortW');
climLong  = strcat(vars,'_LongW');
varsAll = [var climShort climLong addVars];

data = ones(1, length(varsAll));

Weights = array2table(data, 'VariableNames', varsAll);

disp('  Saving generic Weights.mat...')
save(fullfile(inputDir,'Weights.mat'), 'Weights', '-v7.3','-nocompression');

toc

end
