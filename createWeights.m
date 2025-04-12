function Weights = createWeights(targetVar,climateVars,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

climShort = strcat(climateVars,'_ShortW');
climLong  = strcat(climateVars,'_LongW');
targetW   = strcat(targetVar,'_W');

varsAll = [targetW climShort climLong];

data = ones(1, length(varsAll));

Weights = array2table(data, 'VariableNames', varsAll);

disp('  Saving generic Weights.mat...')
save(fullfile(inputDir,'Weights.mat'), 'Weights', '-v7.3','-nocompression');

end
