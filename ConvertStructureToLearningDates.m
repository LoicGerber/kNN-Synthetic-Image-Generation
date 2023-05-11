function learningDates = ConvertStructureToLearningDates(var,LdateStart,LdateEnd,rawData,climateData,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

tic

% Learning dates - variable to be generated
disp('  Processing ' + var + ' for learningDates...')
learningDates = rawData.(lower(var)+'Index');
% Load data set
targetVarData = rawData.(lower(var));
% extract possible learning dates
[r,~] = find(learningDates>=LdateStart & learningDates<=LdateEnd);
learningDates = learningDates(r);
targetVarData = targetVarData(r);

datesAll = rawData.(lower(var)+'Index');

% Find the last common elements of availabletarget variable dates and learningDates
[~, ~, indexLearningDates] = intersect(datesAll, learningDates);
learningDates = learningDates(indexLearningDates); % KEEP ONLY DATES MATCHING CLIMATE DATA
targetVarData = targetVarData(indexLearningDates);

% Find the last common elements of climate dates and learningDates
[~, ~, indexLearningDates] = intersect(climateData.date, learningDates);
learningDates = learningDates(indexLearningDates); % KEEP ONLY DATES MATCHING CLIMATE DATA
targetVarData = targetVarData(indexLearningDates);

learningDates_table = table(learningDates,targetVarData,'VariableNames',['date',var]);
learningDates = learningDates_table;
disp('  Saving Learning dates, may take a while depending on input size...')
save(fullfile(inputDir,'learningDates.mat'), 'learningDates', '-v7.3','-nocompression');

toc

end
