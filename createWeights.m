function Weights = createWeights(climateVars,metricKNN,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

climShort = strcat(climateVars,'_ShortW');
climLong  = strcat(climateVars,'_LongW');

if metricKNN == 5
    spemHel = ["SPEM" "Hellinger"];
else
    spemHel = [];
end

varsAll = [climShort climLong spemHel];

data = ones(1, length(varsAll));

Weights = array2table(data, 'VariableNames', varsAll);

disp('  Saving generic Weights.mat...')
save(fullfile(inputDir,'Weights.mat'), 'Weights', '-v7.3','-nocompression');

end
