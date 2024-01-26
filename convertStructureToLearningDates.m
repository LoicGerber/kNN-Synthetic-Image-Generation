function learningDates = convertStructureToLearningDates(targetVar,LdateStart,LdateEnd,QdateStart,QdateEnd,rawData,climateData,optimPrep,inputDir)

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
    disp("  Processing '" + targetVarL(i) + "' for learningDates...")
    learningDates = rawData.(lower(targetVarL(i))+'Index');
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
for i = 1:numel(targetVarL)
    data = rawData.(lower(targetVarL(i))+'Index');
    [r,~] = find(data>=commonDates(1) & data<=commonDates(end));
    % Load data set
    targetVarData = rawData.(lower(targetVarL(i)));
    targetVarDataAll = [targetVarDataAll targetVarData(r)];
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
learningDatesTable = table(learningDates,targetVarDataAll,'VariableNames',["date" targetVarL']);
% if numel(targetVar) == 1
%     learningDates_table = table(learningDates,targetVarDataAll,'VariableNames',['date',targetVar]);
% elseif numel(targetVar) == 2
%     learningDates_table = table(learningDates,targetVarDataAll(:,1),targetVarDataAll(:,2),'VariableNames',['date',targetVar(1),targetVar(2)]);
% else
%     error('Change ConvertStructureToLearningDates to allow more variable names')
% end
learningDates = learningDatesTable;
disp('  Saving Learning dates, may take a while depending on input size...')
save(fullfile(inputDir,'learningDates.mat'), 'learningDates', '-v7.3','-nocompression');

end
