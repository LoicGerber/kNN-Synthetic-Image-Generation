function [queryDates,learningDates,refValidation] = convertStructureToQueryDates(targetVar,QdateStart,QdateEnd,learningDates,climateData,maxThreshold,validationPrep,optimPrep,outputTime,inputDir,outputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% Query dates - variable to be generated
disp('  Processing queryDates for all target variables...')

targetVarL = lower(targetVar);
datesAll = climateData.date;
learningDatesDate = learningDates.date;

%imgLength = size(learningDates{1,2}{1,1},1);
%imgWidth  = size(learningDates{1,2}{1,1},2);

% Query dates and adapt Learning dates if Validation ON
if validationPrep == false && optimPrep == false % VALIDATION OFF
    % Query dates are all dates in query window, without dates in Learning dates
    if outputTime == 1 % daily
        % Select the dates that are not in learningDates
        r = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        if min(datesAll)>QdateStart
            error('Climate data first date > Query period start')
        elseif max(datesAll)<QdateEnd
            error('Climate data last date < Query period end')
        end
        queryDates = datesAll(r);
        queryDates = setdiff(queryDates, learningDatesDate);
        if isempty(queryDates)
            error('Query dates match with learning dates, nothing to generate')
        end
    elseif outputTime == 2 % monthly
        r = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        if min(datesAll)>QdateStart
            error('Climate data first date > Query period start')
        elseif max(datesAll)<QdateEnd
            error('Climate data last date < Query period end')
        end
        queryDates = datesAll(r);
        % Convert dailyDates to a matrix of year, month, and day components
        dateVec    = datevec(datetime(datesAll,'ConvertFrom','yyyyMMdd'));
        % Find the indices of the dates where the day component is the last day of the month
        lastDays   = find(dateVec(:,3) == eomday(dateVec(:,1), dateVec(:,2)));
        % Select the dates that are not in learningDates
        queryDates = setdiff(queryDates(lastDays), learningDatesDate);
        if isempty(queryDates)
            error('Query dates match with learning dates, nothing to generate')
        end
    elseif outputTime == 3 % dekadal
        r = find(datesAll >= QdateStart & datesAll <= QdateEnd);
        if min(datesAll) > QdateStart
            error('Climate data first date > Query period start')
        elseif max(datesAll) < QdateEnd
            error('Climate data last date < Query period end')
        end
        queryDates = datesAll(r);
        % Convert dailyDates to a matrix of year, month, and day components
        dateVec = datevec(datetime(queryDates, 'ConvertFrom', 'yyyyMMdd'));
        % Find the indices of the dates for the 10th, 20th, and the last day of the month
        selectedDays = [];
        selectedDays = [selectedDays; find(dateVec(:,3) == 10)];
        selectedDays = [selectedDays; find(dateVec(:,3) == 20)];
        selectedDays = [selectedDays; find(dateVec(:,3) == eomday(dateVec(:, 1), dateVec(:, 2)))];
        selectedDays = sort(selectedDays);
        % Select the dates that are not in learningDates
        queryDates = setdiff(queryDates(selectedDays), learningDatesDate);
        if isempty(queryDates)
            error('Query dates match with learning dates, nothing to generate...')
        end
    else
        error('Invalid outputTime value')
    end
    refValidation = [];
elseif validationPrep == true || optimPrep == true % validation or optimPrep ON
    % Query dates are all dates in query window, replacing dates in Learning dates
    for j = 1:numel(targetVarL)
        if ~exist(outputDir,'dir')
            mkdir(outputDir)
        end
    end
    if outputTime == 1 % daily
        % Select the dates that are not in learningDates
        r = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        if min(datesAll)>QdateStart
            error('Climate data first date > Query period start')
        elseif max(datesAll)<QdateEnd
            error('Climate data last date < Query period end')
        end
        queryDates = datesAll(r);
    elseif outputTime == 2 % monthly
        r = find(datesAll>=QdateStart & datesAll<=QdateEnd);
        if min(datesAll)>QdateStart
            error('Climate data first date > Query period start')
        elseif max(datesAll)<QdateEnd
            error('Climate data last date < Query period end')
        end
        queryDates = datesAll(r);
        % Convert dailyDates to a matrix of year, month, and day components
        dateVec    = datevec(datetime(queryDates,'ConvertFrom','yyyyMMdd'));
        % Find the indices of the dates where the day component is the last day of the month
        lastDays   = find(dateVec(:,3) == eomday(dateVec(:,1), dateVec(:,2)));
        % Select the dates that are not in learningDates
        queryDates = queryDates(lastDays);
    elseif outputTime == 3 % dekadal
        r = find(datesAll >= QdateStart & datesAll <= QdateEnd);
        if min(datesAll) > QdateStart
            error('Climate data first date > Query period start')
        elseif max(datesAll) < QdateEnd
            error('Climate data last date < Query period end')
        end
        queryDates = datesAll(r);
        % Convert dailyDates to a matrix of year, month, and day components
        dateVec = datevec(datetime(queryDates, 'ConvertFrom', 'yyyyMMdd'));
        % Find the indices of the dates for the 10th, 20th, and the last day of the month
        selectedDays = [];
        selectedDays = [selectedDays; find(dateVec(:,3) == 10)];
        selectedDays = [selectedDays; find(dateVec(:,3) == 20)];
        selectedDays = [selectedDays; find(dateVec(:,3) == eomday(dateVec(:, 1), dateVec(:, 2)))];
        selectedDays = sort(selectedDays);
        % Select the dates that are not in learningDates
        queryDates = queryDates(selectedDays);
    else
        error('Invalid outputTime value')
    end
    % Define learningDates as itself minus the query dates
    ismem = ismember(learningDatesDate, queryDates);
    learningDataValidation  = learningDates(~ismem,:);
    learningDatesValidation = learningDatesDate(ismem);
    referenceValidation = {};
    for j = 1:numel(targetVarL)
        referenceValidation = [referenceValidation table2cell(learningDates(ismem,targetVarL(j)))];
        imagesRefValidation = nan(size(referenceValidation{1,j},1),size(referenceValidation{1,j},2),size(learningDatesValidation,1));
        % Create matrix of reference dates
        for i = 1:size(learningDatesValidation,1)
            imagesRefValidation(:,:,i) = referenceValidation{i,j};
        end
        refValidation.(targetVar(j)) = single(imagesRefValidation);
    end
    refValidation.date = learningDatesValidation;
    disp('  Saving refValidation.mat file...')
    save(fullfile(inputDir,'refValidation.mat'), 'refValidation', '-v7.3','-nocompression');
    learningDates = learningDataValidation;
end

% Select closest targetVar index for each Query date
nearestIdx = nan(size(queryDates));
if validationPrep == false && optimPrep == false % validation OFF
    for i = 1:numel(queryDates)
        %[nearest, nearestIdx(i)] = min(abs(learningDatesDate - queryDates(i)));  % find index of closest date
        [nearest, nearestIdx(i)] = min(abs(datetime(learningDatesDate,'ConvertFrom','yyyymmdd') - datetime(queryDates(i),'ConvertFrom','yyyyMMdd')));
        if days(nearest) > maxThreshold %%% MAX THRESHOLD <--------------------------------------------------------------------------------------------------------------------------------
            nearestIdx(i) = nan;
        end
    end
    % Match the dates
    matchedTargetVarDates = [queryDates, nan(size(queryDates))];
    for i = 1:length(nearestIdx)
        if ~isnan(nearestIdx(i))
            matchedTargetVarDates(i,2) = learningDatesDate(nearestIdx(i));
        end
    end
    % Assign closest targetVar map to each Query date
    %     try
    %         matchedTargetVarTable = table('Size',[size(matchedTargetVarDates,1),numel(targetVar)+1], 'VariableTypes',{'double', 'cell'});
    %     catch
    %         try
    %             matchedTargetVarTable = table('Size',[size(matchedTargetVarDates,1),numel(targetVar)+1], 'VariableTypes',{'double', 'cell', 'cell'});
    %         catch
    %             error('Adapt ConvertStructureToQueryDates function to allow more variables')
    %         end
    %     end
    for j = 1:numel(targetVarL)
        targetVarData = learningDates.(targetVarL(j));
        % Loop through the matched dates
        for i = 1:size(matchedTargetVarDates, 1)
            % Get the date to match
            matchDate = matchedTargetVarDates(i, 2);
            % If a match was found, add the date and data to the output table
            if j == 1, matchedTargetVarTable{i, 1}   = queryDates(i); end
            if ~isnan(matchDate)
                matchedTargetVarTable{i, j+1} = targetVarData(nearestIdx(i));
            else
                matchedTargetVarTable{i, j+1} = {nan};
            end
        end
    end
elseif validationPrep == true || optimPrep == true % validation or optimPrep ON
    matchedTargetVarDates = [queryDates, nan(size(queryDates))];
    %     try
    %         matchedTargetVarTable = table('Size',[size(matchedTargetVarDates,1),numel(targetVar)+1], 'VariableTypes',{'double', 'cell'});
    %     catch
    %         try
    %             matchedTargetVarTable = table('Size',[size(matchedTargetVarDates,1),numel(targetVar)+1], 'VariableTypes',{'double', 'cell', 'cell'});
    %         catch
    %             error('Adapt ConvertStructureToQueryDates function to allow more variables')
    %         end
    %     end
    for j = 1:numel(targetVarL)
        targetVarData = learningDates.(targetVarL(j));
        % Loop through the matched dates
        for i = 1:size(matchedTargetVarDates, 1)
            % fill the table with NaNs the size of the variable to be generated
            if j == 1, matchedTargetVarTable{i, j} = queryDates(i); end
            matchedTargetVarTable{i, j+1} = {nan};
            %matchedTargetVarTable{i, j+1} = {nan(size(targetVarData{1,1}))};
        end
    end
end
% Rename the columns
% try
%     matchedTargetVarTable.Properties.VariableNames = {'Date', convertStringsToChars(targetVar)};
% catch
%     try
%         matchedTargetVarTable.Properties.VariableNames = {'Date', convertStringsToChars(targetVar(1)),convertStringsToChars(targetVar(2))};
%     catch
%         error('Adapt ConvertStructureToQueryDates function to allow more variables')
%     end
% end
matchedTargetVarTable = cell2table(matchedTargetVarTable,"VariableNames",["date" targetVarL']);
queryDates = matchedTargetVarTable;
disp('  Saving Query dates, may take a while depending on input size...')
save(fullfile(inputDir,'queryDates.mat'), 'queryDates', '-v7.3','-nocompression');

end
