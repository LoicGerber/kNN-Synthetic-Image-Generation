function climateData = extractClimateData(vars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,longWindow,validation,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

tic

% Get fieldnames matching the string
matchingFields     = string(fieldnames(rawData));
matchingDataFields = matchingFields(ismember(lower(matchingFields), lower(vars)));
datesFields        = matchingFields(ismember(lower(matchingFields), lower(vars+'index')));

% Create cell array with each variable in a separate column
%num_vars = length(matchingDataFields);
%climateData = cell(numel(output_data.(matchingDataFields(1))), num_vars);
climateData = table();
for i = 1:length(matchingDataFields)
    climateData(:,i) = rawData.(matchingDataFields(i));
end
climateData.Properties.VariableNames = matchingDataFields;
climateData.("date") = rawData.(datesFields(1));
climateData = movevars(climateData,'date','Before',1);

if validation == 1 % validation OFF
    % Take all available climate dates
    [r,~] = find(datetime(climateData.date,'ConvertFrom','yyyymmdd')>=datetime(min(QdateStart,LdateStart),'ConvertFrom','yyyyMMdd')-days(longWindow) ...
        & datetime(climateData.date,'ConvertFrom','yyyymmdd')<=datetime(max(QdateEnd,LdateEnd),'ConvertFrom','yyyymmdd'));
elseif validation == 0 % validation ON
    % Take dates only within Learning dates range
    [r,~] = find(datetime(climateData.date,'ConvertFrom','yyyymmdd')>=datetime(LdateStart,'ConvertFrom','yyyymmdd')-days(longWindow) ...
        & datetime(climateData.date,'ConvertFrom','yyyymmdd')<=datetime(LdateEnd,'ConvertFrom','yyyymmdd'));
else
    error('Validation value is either 0 (OFF) or 1 (ON)...')
end
climateData = climateData(r,:);

disp('Saving climateData table...')
allVarsSave = fullfile(inputDir, 'climateData.mat');
save(allVarsSave, 'climateData', '-v7.3','-nocompression');

toc

end
