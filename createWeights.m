function Weights = createWeights(targetVar,climateVars,addVars,inputDir)

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
if ~isempty(addVars)
    addVarsW  = strcat(addVars,'_W');
else
    addVarsW = [];
end
varsAll = [targetW climShort climLong addVarsW];

data = ones(1, length(varsAll));

Weights = array2table(data, 'VariableNames', varsAll);

disp('  Saving generic Weights.mat...')
save(fullfile(inputDir,'Weights.mat'), 'Weights', '-v7.3','-nocompression');

end
