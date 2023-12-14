function additionalVars = extractAdditionalVars(addVars,rawData,climateData,QdateStart,QdateEnd,LdateStart,LdateEnd,maxThreshold,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

if ~isempty(addVars)
    % Get fieldnames matching the string
    matchingFields     = string(fieldnames(rawData));
    matchingDataFields = matchingFields(ismember(lower(matchingFields), lower(addVars)));
    datesFields        = matchingFields(ismember(lower(matchingFields), lower(addVars+'index')));

    allDates           = climateData.date;

    additionalVars = table();
    datesAddVars = [];
    for i = 1:length(matchingDataFields)
        datesAddVars = [datesAddVars; rawData.(datesFields(i))];
    end
    datesAddVars = sort(unique(datesAddVars));
    additionalVars.("date") = datesAddVars;

    for i = 1:length(matchingDataFields)
        additionalVars.(matchingDataFields(i)) = cell(length(datesAddVars),1);
        datesIdx = ismember(datesAddVars,rawData.(datesFields(i)));
        additionalVars.(matchingDataFields(i))(datesIdx) = rawData.(matchingDataFields(i));
        %additionalVars.(matchingDataFields(i))(~datesIdx) = {single(nan(size(rawData.(matchingDataFields(i)){1,1})))};
    end

    % Take all available climate dates
    r = find(datetime(additionalVars.date,'ConvertFrom','yyyymmdd')>=datetime(min(QdateStart,LdateStart),'ConvertFrom','yyyyMMdd') ...
        & datetime(additionalVars.date,'ConvertFrom','yyyymmdd')<=datetime(max(QdateEnd,LdateEnd),'ConvertFrom','yyyymmdd'));
    if min(datetime(additionalVars.date,'ConvertFrom','yyyymmdd'))>datetime(LdateStart,'ConvertFrom','yyyymmdd')
        warning('Additional data first date > Learning period start')
    end
    if max(datetime(additionalVars.date,'ConvertFrom','yyyymmdd'))<datetime(LdateEnd,'ConvertFrom','yyyymmdd')
        warning('Additional data last date < Learning period end')
    end
    if ~(min(datetime(additionalVars.date,'ConvertFrom','yyyymmdd'))>datetime(LdateEnd,'ConvertFrom','yyyymmdd') || ...
            max(datetime(additionalVars.date,'ConvertFrom','yyyymmdd'))<datetime(LdateStart,'ConvertFrom','yyyymmdd'))
        additionalVars = additionalVars(r,:);
        datesUnique = cell(1,numel(matchingDataFields));
        for j = 1:numel(matchingDataFields)
            rU = ismember(rawData.(datesFields(j)),additionalVars.date);
            datesUnique{1,j} = rawData.(datesFields(j))(rU);
            dataUnique{1,j} = rawData.(matchingDataFields(j))(rU);
        end

        % Select closest addVar index for each Query date
        nearestIdx = nan(numel(allDates),numel(matchingDataFields));
        for j = 1:numel(matchingDataFields)
            for i = 1:numel(allDates)
                %[nearest, nearestIdx(i)] = min(abs(learningDatesDate - queryDates(i)));  % find index of closest date
                [nearest, nearestIdx(i,j)] = min(abs(datetime(cell2mat(datesUnique(j)),'ConvertFrom','yyyymmdd') - datetime(allDates(i),'ConvertFrom','yyyymmdd')));
                if days(nearest) > maxThreshold
                    nearestIdx(i,j) = nan;
                end
            end
        end
        % Match the dates
        matchedAddVarDates = allDates;
        for j = 1:numel(matchingDataFields)
            datesU = cell2mat(datesUnique(j));
            matchedAddVarDates = [matchedAddVarDates, nan(size(allDates))];
            for i = 1:length(nearestIdx)
                if ~isnan(nearestIdx(i,j))
                    matchedAddVarDates(i,j+1) = datesU(nearestIdx(i,j));
                end
            end
        end
        % Assign closest addVar map to each Query date
        %     try
        %         matchedAddVarTable = table('Size',[size(matchedAddVarDates,1),numel(matchingDataFields)+1], 'VariableTypes',{'double', 'cell'});
        %     catch
        %         try
        %             matchedAddVarTable = table('Size',[size(matchedAddVarDates,1),numel(matchingDataFields)+1], 'VariableTypes',{'double', 'cell', 'cell'});
        %         catch
        %             error('Adapt extractAdditionalVars function to allow more variables')
        %         end
        %     end
        for j = 1:numel(matchingDataFields)
            addVarData = additionalVars.(matchingDataFields(j));
            dataU = dataUnique{j};
            % Loop through the matched dates
            for i = 1:size(matchedAddVarDates, 1)
                % Get the date to match
                matchDate = matchedAddVarDates(i, j+1);
                % If a match was found, add the date and data to the output table
                if j == 1, matchedAddVarTable{i, 1} = allDates(i); end
                if ~isnan(matchDate)
                    %matchedAddVarTable{i, j+1} = addVarData(nearestIdx(i,j));
                    matchedAddVarTable{i, j+1} = dataU{nearestIdx(i,j)};
                    %else
                    %matchedAddVarTable{i, j+1} = {nan(size(addVarData{1,1}))};
                end
            end
        end
        % Rename the columns
        %     try
        %         matchedAddVarTable.Properties.VariableNames = {'date', convertStringsToChars(matchingDataFields)};
        %     catch
        %         try
        %             matchedAddVarTable.Properties.VariableNames = {'date', convertStringsToChars(matchingDataFields(1)),convertStringsToChars(matchingDataFields(2))};
        %         catch
        %             error('Adapt extractAdditionalVars function to allow more variables')
        %         end
        %     end
        matchedAddVarTable = cell2table(matchedAddVarTable,"VariableNames",["date" matchingDataFields']);
        additionalVars = matchedAddVarTable;
    else
        warning('Additional data dates not within Learning days range...')
        additionalVars = [];
    end
else
    additionalVars = [];
end

disp('  Saving additionalVars table...')
allVarsSave = fullfile(inputDir, 'additionalVars.mat');
save(allVarsSave, 'additionalVars', '-v7.3','-nocompression');


end

