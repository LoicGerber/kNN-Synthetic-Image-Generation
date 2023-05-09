function additionalVars = extractAdditionalVars(addVars,rawData,QdateStart,QdateEnd,LdateStart,LdateEnd,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

tic

if ~isempty(addVars)
    % Get fieldnames matching the string
    matchingFields     = string(fieldnames(rawData));
    matchingDataFields = matchingFields(ismember(lower(matchingFields), lower(addVars)));
    datesFields        = matchingFields(ismember(lower(matchingFields), lower(addVars+'index')));

    additionalVars = table();
    for i = 1:length(matchingDataFields)
        additionalVars(:,i) = rawData.(matchingDataFields(i));
    end
    additionalVars.Properties.VariableNames = matchingDataFields;
    additionalVars.("date") = rawData.(datesFields(1));
    additionalVars = movevars(additionalVars,'date','Before',1);

    % Take all available climate dates
    [r,~] = find(datetime(additionalVars.date,'ConvertFrom','yyyymmdd')>=datetime(min(QdateStart,LdateStart),'ConvertFrom','yyyyMMdd') ...
        & datetime(additionalVars.date,'ConvertFrom','yyyymmdd')<=datetime(max(QdateEnd,LdateEnd),'ConvertFrom','yyyymmdd'));
    additionalVars = additionalVars(r,:);
else
    additionalVars = [];
end

disp('Saving additionalVars table...')
allVarsSave = fullfile(inputDir, 'additionalVars.mat');
save(allVarsSave, 'additionalVars', '-v7.3','-nocompression');
disp('additionalVars table saved')

toc

end

