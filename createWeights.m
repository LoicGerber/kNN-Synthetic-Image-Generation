function Weights = createWeights(climateVars,metricKNN,useDOY,inputDir)

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
if useDOY
    doy = "DOY";
else
    doy = [];
end

varsAll = [climShort climLong doy spemHel];

data = ones(1, length(varsAll));

Weights = array2table(data, 'VariableNames', varsAll);

disp('  Saving generic Weights.mat...')
save(fullfile(inputDir,'Weights.mat'), 'Weights', '-v7.3','-nocompression');

end
