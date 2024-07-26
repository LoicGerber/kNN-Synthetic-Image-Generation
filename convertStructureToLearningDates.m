function learningDates = convertStructureToLearningDates(targetVar,LdateStart,LdateEnd,QdateStart,QdateEnd,rawData,climateData,targetDim,optimPrep,inputDir,saveMats)

%
%
%
% REDO DOCUMENTATION
%
%
%

targetVarL = lower(targetVar);
commonDates = [];
for i = 1:numel(targetVarL)
    % Learning dates - variable to be generated
    disp("  Processing '" + targetVar(i) + "' for learningDates...")
    learningDates = rawData.(strcat(lower(targetVarL(i)),'Index'));
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
if targetDim ~= 1
    targetVarDataAll = {};
    for i = 1:numel(targetVarL)
        data = rawData.(strcat(lower(targetVarL(i)),'Index'));
        [r,~] = find(ismember(data,commonDates));
        % Load data set
        targetVarData = rawData.(lower(targetVarL(i)));
        targetVarDataAll = [targetVarDataAll targetVarData(r)];
    end
else
    for i = 1:numel(targetVarL)
        data = rawData.(strcat(lower(targetVarL(i)),'Index'));
        [r,~] = find(ismember(data,commonDates));
        % Load data set
        targetVarData = rawData.(lower(targetVarL(i)));
        targetVarData = targetVarData(r);
    end
end
% Keep only the target variable data corresponding to the common dates
r = find(learningDates>=LdateStart & learningDates<=LdateEnd);
if optimPrep == true && (LdateStart ~= QdateStart)
    % if in preparation for optimisation, include Query dates in Learning dates
    rQ = find(learningDates>=QdateStart & learningDates<=QdateEnd);
    r = unique([r; rQ]);
end
if min(learningDates)>LdateStart
    warning(['Target data first date > Learning period start (',num2str(min(learningDates)),' vs ',num2str(LdateStart),')'])
elseif max(learningDates)<LdateEnd
    warning(['Target data last date < Learning period end (',num2str(max(learningDates)),' vs ',num2str(LdateEnd),')'])
end
learningDates = learningDates(r);
if targetDim ~= 1
    targetVarData = targetVarDataAll(r);
else
    targetVarData = targetVarData(r);
end
datesAll = commonDates;

% Find the last common elements of target variable dates and learningDates
[~, ~, indexLearningDates] = intersect(datesAll, learningDates);
learningDates = learningDates(indexLearningDates); % KEEP ONLY DATES MATCHING CLIMATE DATA
if targetDim ~= 1
    targetVarDataAll = targetVarDataAll(indexLearningDates,:);
else
    targetVarDataAll = targetVarData(indexLearningDates);
end

% Find the last common elements of climate dates and learningDates
[~, ~, indexLearningDates] = intersect(climateData.date, learningDates);
learningDates = learningDates(indexLearningDates); % KEEP ONLY DATES MATCHING CLIMATE DATA
targetVarDataAll = targetVarDataAll(indexLearningDates,:);
learningDatesTable = table(learningDates,'VariableNames',"date");
for i = 1:numel(targetVar)
    learningDatesTable.(targetVarL(i)) = targetVarDataAll(:,i);
end
learningDates = learningDatesTable;

if saveMats == true
    disp('  Saving Learning dates, may take a while depending on input size...')
    save(fullfile(inputDir,'learningDates.mat'), 'learningDates', '-v7.3','-nocompression');
end

end
