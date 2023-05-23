function allWeights = checkWeights(Weights,learningDates)

%
%
%
% REDO DOCUMENTATION
%
%
%

varNames       = learningDates.Properties.VariableNames;       % Extract columns name
varNames       = varNames(~ismember(varNames,{'Dates','Date'})); % Remove the 'Dates' variable name
uniqueVarNames = cellfun(@(x) strsplit(x, '-'), varNames, 'UniformOutput', false);
uniqueVarNames = unique(cellfun(@(x) x{1}, uniqueVarNames, 'UniformOutput', false));

% If Weights Table is inputed
if ~isempty(Weights)
    disp('Weights matrix inputed!')
    weightsFields = Weights.Properties.VariableNames;
    % Loop over the fields
    for i = 1:numel(weightsFields)
        % Get the name of the current field
        currentField = weightsFields{i};
        % Get the data associated with the current field
        currentFieldData = Weights.(currentField);
        % Create a new variable with the name of the current field and populate it with the associated data
        eval([currentField, '=currentFieldData;']);
    end
else % if table is not inputed, creates weights for all variables contained in the input file
    disp('No Weights matrix inputed, creating generic weights equal to 1 for all variables...')
    weightsFields = uniqueVarNames;
    Weights = struct();
    i = 1;
    while i <= length(weightsFields)
        if shortWindow ~= 0
            if ismember(weightsFields{i}, vars)
                newVarName = weightsFields{i};
                weightsFields = [weightsFields(1:i-1) {newVarName} weightsFields(i:end)];
                % Modify the duplicated name with 'wClose' prefix
                weightsFields{i+1} = ['Close', newVarName];
            end
        end
        % Add 'w' prefix to original variable name
        weightsFields{i} = ['w', weightsFields{i}];
        % Increment i by 1 to move to the next variable name
        i = i + 1;
    end
    for j = 1:numel(weightsFields)
        % Get the name of the current variable
        currentField = weightsFields{j};
        % Get the data associated with the current field
        Weights.(currentField) = 1;
    end
end

% CHECK THAT ALL WEIGHTS MATCH VARIABLES! IF NOT, DISPLAY ERROR MESSAGE IF WEIGHT IS MISSING
notContained = {};
for i = 1:length(uniqueVarNames)
    if ~any(contains(weightsFields, uniqueVarNames{i}))
        notContained{end+1} = uniqueVarNames{i};
    end
end
if ~isempty(notContained)
    error(['The following variables do not have associated weights in weightsFields: ', strjoin(notContained, ', ')]);
end

disp('Assigning weights to AllWeights matrix...')
% Makes a weights matrix
allWeights = zeros(size(learning_dates,1),size(learning_dates,2));
% Create an empty structure
variableIndices = struct();
% Loop over each unique variable name
for i = 1:length(uniqueVarNames)
    % Find the indices of columns containing the current variable name
    currentVarName = uniqueVarNames{i};
    currentIndices = find(contains(varNames, uniqueVarNames{i}));
    if shortWindow ~= 0
        if strcmpi(currentVarName, var) % For variable to be generated
            % Assign the current indices to the corresponding field in the structure
            variableIndices.(uniqueVarNames{i}) = currentIndices;
        else % For climate variables
            % Assign the current indices to the corresponding field in the structure
            variableIndices.(uniqueVarNames{i}) = currentIndices(1:end-shortWindow);
            variableIndices.(strcat('Close',uniqueVarNames{i})) = currentIndices((end-shortWindow)+1:end);
        end
    else
        % Assign the current indices to the corresponding field in the structure
        variableIndices.(uniqueVarNames{i}) = currentIndices(1:end);
    end
end

% Asign weights to correct indices of the weight matrix
variableIdNames = fieldnames(variableIndices)';
for i = 1:length(variableIdNames)
    % Find the indices of columns containing the current variable name
    currentIndices = variableIndices.(variableIdNames{i});
    % Get the corresponding weights for the current variable
    currentWeights = Weights.(strcat('w', variableIdNames{i}));
    % Assign the weights to the corresponding columns in AllWeights
    allWeights(:,currentIndices) = currentWeights.';
end

SumWeights = sum(allWeights(1,:));
allWeights = allWeights./SumWeights; % normalisation of allWeights

end
