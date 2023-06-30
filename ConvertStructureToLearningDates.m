function learningDates = ConvertStructureToLearningDates(var,LdateStart,LdateEnd,QdateStart,QdateEnd,rawData,climateData,optimPrep,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

commonDates = [];
for i = 1:numel(var)
    % Learning dates - variable to be generated
    disp("  Processing '" + var(i) + "' for learningDates...")
    learningDates = rawData.(lower(var(i))+'Index');
    % Extract the common dates
    if isempty(commonDates)
        % For the first variable, set the common dates as all available learning dates
        commonDates = learningDates;
    else
        % For subsequent variables, keep only the dates that exist in commonDates and learningDates
        commonDates = intersect(commonDates, learningDates);
    end
end
learningDates = commonDates;
targetVarDataAll = {};
for i = 1:numel(var)
    data = rawData.(lower(var(i))+'Index');
    [r,~] = find(data>=commonDates(1) & data<=commonDates(end));
    % Load data set
    targetVarData = rawData.(lower(var(i)));
    targetVarDataAll = [targetVarDataAll targetVarData(r)];
end
% Keep only the target variable data corresponding to the common dates
[r,~] = find(learningDates>=LdateStart & learningDates<=LdateEnd);
if optimPrep == true && (LdateStart ~= QdateStart)
    % if in preparation for optimisation, include Query dates in Learning dates
    [rQ,~] = find(learningDates>=QdateStart & learningDates<=QdateEnd);
    r = unique([r; rQ]);
end
learningDates = learningDates(r);
targetVarData = targetVarDataAll(r);
    
datesAll = commonDates;

% Find the last common elements of target variable dates and learningDates
[~, ~, indexLearningDates] = intersect(datesAll, learningDates);
learningDates = learningDates(indexLearningDates); % KEEP ONLY DATES MATCHING CLIMATE DATA
targetVarDataAll = targetVarDataAll(indexLearningDates,:);

% Find the last common elements of climate dates and learningDates
[~, ~, indexLearningDates] = intersect(climateData.date, learningDates);
learningDates = learningDates(indexLearningDates); % KEEP ONLY DATES MATCHING CLIMATE DATA
targetVarDataAll = targetVarDataAll(indexLearningDates,:);

if numel(var) == 1
    learningDates_table = table(learningDates,targetVarDataAll,'VariableNames',['date',var]);
elseif numel(var) == 2
    learningDates_table = table(learningDates,targetVarDataAll(:,1),targetVarDataAll(:,2),'VariableNames',['date',var(1),var(2)]);
else
    error('Change ConvertStructureToLearningDates to allow more variable names')
end
learningDates = learningDates_table;
disp('  Saving Learning dates, may take a while depending on input size...')
save(fullfile(inputDir,'learningDates.mat'), 'learningDates', '-v7.3','-nocompression');

end
