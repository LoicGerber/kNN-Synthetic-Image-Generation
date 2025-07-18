function climateData = extractClimateData(climateVars,rawData,normMethods,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,inputDir,saveMats)

%
%
%
% REDO DOCUMENTATION
%
%
%

% Get fieldnames matching the string
matchingFields     = string(fieldnames(rawData));
matchingDataFields = matchingFields(ismember(lower(matchingFields), lower(climateVars)));
datesFields        = matchingFields(ismember(lower(matchingFields), lower(climateVars+'index')));

% Create cell array with each variable in a separate column
%num_vars = length(matchingDataFields);
%climateData = cell(numel(output_data.(matchingDataFields(1))), num_vars);
climateData = table();
for i = 1:length(matchingDataFields)
    currentName = matchingDataFields(i);
    currentIdx  = find(strcmpi(climateVars,currentName));
    if normMethods(currentIdx) == 0 %No normalisation
        climateData(:,i) = rawData.(matchingDataFields(i));
    elseif normMethods(currentIdx) == 1 %MinMax
        % Flatten the cell array into a single numeric array
        numMatrices = numel(rawData.(matchingDataFields(i)));
        flattenedData = nan(numMatrices, numel(rawData.(matchingDataFields(i)){1})); % Initialize a matrix to store flattened data
        for j = 1:numMatrices
            flattenedData(j, :) = rawData.(matchingDataFields(i)){j}(:); % Flatten each matrix and store in the matrix
        end
        % Find the minimum and maximum values of the entire numeric array
        minValue = min(flattenedData(:));
        maxValue = max(flattenedData(:));
        % Normalize each matrix in the cell array
        normalizedCellArray = cell(numMatrices, 1);
        for j = 1:numMatrices
            normalizedCellArray{j} = (rawData.(matchingDataFields(i)){j} - minValue) / (maxValue - minValue);
        end
        % climateData is normalised
        climateData(:,i) = normalizedCellArray;
    elseif normMethods(currentIdx) == 2 %Q5-Q95
        % Flatten the cell array into a single numeric array
        numMatrices = numel(rawData.(matchingDataFields(i)));
        flattenedData = nan(numMatrices, numel(rawData.(matchingDataFields(i)){1})); % Initialize a matrix to store flattened data
        for j = 1:numMatrices
            flattenedData(j, :) = rawData.(matchingDataFields(i)){j}(:); % Flatten each matrix and store in the matrix
        end
        % Find the Q5 and Q95 values of the entire numeric array
        Q5  = quantile(flattenedData(:),0.05);
        Q95 = quantile(flattenedData(:),0.95);
        % Normalize each matrix in the cell array
        normalizedCellArray = cell(numMatrices, 1);
        for j = 1:numMatrices
            normRawData = rawData.(matchingDataFields(i)){j} > Q95;
            rawData.(matchingDataFields(i)){j}(normRawData) = Q95;
            normRawData = rawData.(matchingDataFields(i)){j} < Q5;
            rawData.(matchingDataFields(i)){j}(normRawData) = Q5;
            normalizedCellArray{j} = (rawData.(matchingDataFields(i)){j} - Q5) / (Q95 - Q5);
        end
        % climateData is normalised
        climateData(:,i) = normalizedCellArray;
    elseif normMethods(currentIdx) == 3 % log
        % Transform data to log(data+1)
        %climateData(:,i) = cellfun(@(x) log(x + 1) .* (x > 0),rawData.(matchingDataFields(i)),'UniformOutput',false);
        climateData(:,i) = cellfun(@(x) log(x + 1), rawData.(matchingDataFields(i)),'UniformOutput',false);
    elseif normMethods(currentIdx) == 4 % 
        % Transform data to log(data+1), only for data > 0
        climateData(:,i) = cellfun(@(x) log(x + 1) .* (x > 0),rawData.(matchingDataFields(i)),'UniformOutput',false);
        %climateData(:,i) = cellfun(@(x) log(x + 1), rawData.(matchingDataFields(i)),'UniformOutput',false);
    end
    % climateData is not normalised
    %climateData(:,i) = rawData.(matchingDataFields(i));
end
climateData.Properties.VariableNames = matchingDataFields;
climateData.("date") = rawData.(datesFields(1));
climateData = movevars(climateData,'date','Before',1);

if length(QdateStart) > 1
    QdateStart = QdateStart(1);
    QdateEnd   = QdateEnd(end);
end

[rLD,~] = find(datetime(climateData.date,'ConvertFrom','yyyymmdd')>=datetime(LdateStart,'ConvertFrom','yyyymmdd')-days(longWindow) ...
    & datetime(climateData.date,'ConvertFrom','yyyymmdd')<=datetime(LdateEnd,'ConvertFrom','yyyymmdd'));
if min(datetime(climateData.date,'ConvertFrom','yyyymmdd'))>datetime(LdateStart,'ConvertFrom','yyyymmdd')
    error('Climate data first date > Learning period start')
elseif max(datetime(climateData.date,'ConvertFrom','yyyymmdd'))<datetime(LdateEnd,'ConvertFrom','yyyymmdd')
    error('Climate data last date < Learning period end')
end

[rQD,~] = find(datetime(climateData.date,'ConvertFrom','yyyymmdd')>=datetime(QdateStart,'ConvertFrom','yyyymmdd')-days(longWindow) ...
    & datetime(climateData.date,'ConvertFrom','yyyymmdd')<=datetime(QdateEnd,'ConvertFrom','yyyymmdd'));
if min(datetime(climateData.date,'ConvertFrom','yyyymmdd'))>datetime(QdateStart,'ConvertFrom','yyyymmdd')
    error('Climate data first date > Query period start')
elseif max(datetime(climateData.date,'ConvertFrom','yyyymmdd'))<datetime(QdateEnd,'ConvertFrom','yyyymmdd')
    error('Climate data last date < Query period end')
end
r = unique([rLD; rQD]);

climateData = climateData(r,:);

if saveMats == true
    disp('  Saving climateData table...')
    allVarsSave = fullfile(inputDir, 'climateData.mat');
    save(allVarsSave, 'climateData', '-v7.3','-nocompression');
end

end
