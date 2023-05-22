function sortedDates = KNNSortingOptim(var,addVars,shortWindow,Weights,nbImages,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% 1 Q date, 2 L dates, 3 target distance, 4 addVars distances, 5 climate distances, 6 std
distances = load(fullfile(inputDir,'KNNDistances.mat'));
distances = distances.sortedDates;

% Assign different weights
idxTarget     = contains(Weights.Properties.VariableNames,var);
weightsTarget = table2array(Weights(:,idxTarget));
idxShort      = contains(Weights.Properties.VariableNames,'Short');
weightsShort  = table2cell(Weights(:,idxShort));
idxLong       = contains(Weights.Properties.VariableNames,'Long');
weightsLong   = table2cell(Weights(:,idxLong));
if ~isempty(addVars)
    idxAddVars     = contains(Weights.Properties.VariableNames,addVars);
    weightsAddVars = table2cell(Weights(:,idxAddVars));
else
    weightsAddVars = [];
end

totQDates   = size(distances,1);
totLDates   = size(distances{1,2},1);

targetDistance  = cell(totQDates, 1);
addVarsDistance = cell(totQDates, 2);
climateDistance = cell(totQDates, 2);
stdDistance     = cell(totQDates, 1);

sortedDates = cell(totQDates, 1);
sortedData  = cell(totQDates, 1);
sortedDist  = cell(totQDates, 1);
sortedStd   = cell(totQDates, 1);

for qd = 1:totQDates
    for ld = 1:totLDates
        currentLDate = distances{1,2}{ld,1};
        % Target variable comparison
        targetDistance{ld,1} = distances{qd,3}{ld,1};
        targetDistance{ld,1} = targetDistance{ld,1}.*weightsTarget;

        % Additional variable comparison
        if ~isempty(addVars)
            if numel(addVars) == 1
                addVarsDistance{ld,1} = distances{qd,4}{ld,1};
                addVarsDistance{ld,1} = addVarsDistance{ld,1} .* cell2mat(weightsAddVars);
                addVarsDistance{ld,2} = addVarsDistance{ld,2} .* cell2mat(weightsAddVars);
            else
                addVarsDistance{:,ld} = distances{:,ld};
                addVarsDistance{:,ld} = num2cell(cell2mat(addVarsDistance{:,ld}) .* cell2mat(weightsAddVars));
                addVarsDistance{ld,2} = num2cell(cell2mat(addVarsDistance{ld,2}) .* cell2mat(weightsAddVars));
            end
        else
            addVarsDistance{ld,1} = 0;
            addVarsDistance{ld,2} = 0;
        end

        % Climate distance
        climateDistance{ld,1} = currentLDate;
        climateDistance{ld,2} = distances{qd,5}{ld,1};
        climateDistance{ld,2}(1:shortWindow,:)     = num2cell(cell2mat(climateDistance{ld,2}(1:shortWindow,:))     .* cell2mat(weightsShort));
        climateDistance{ld,2}(shortWindow+1:end,:) = num2cell(cell2mat(climateDistance{ld,2}(shortWindow+1:end,:)) .* cell2mat(weightsLong));
        climateDistance{ld,2} = sum(cellfun(@double,climateDistance{ld,2}),1,'omitnan');
        climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+targetDistance{ld,1}+addVarsDistance{ld,1};
        
        % Std climate distance
        stdDistance{ld,1} = distances{qd,6}{ld,1};
        stdDistance{ld,1}(1:shortWindow,:)     = num2cell(cell2mat(stdDistance{ld,1}(1:shortWindow,:))     .* cell2mat(weightsShort));
        stdDistance{ld,1}(shortWindow+1:end,:) = num2cell(cell2mat(stdDistance{ld,1}(shortWindow+1:end,:)) .* cell2mat(weightsLong));
        stdDistance{ld,1} = sum(cellfun(@double,stdDistance{ld,1}),1,'omitnan');
        stdDistance{ld,1} = sum(stdDistance{ld,1},2,'omitnan')+addVarsDistance{ld,2};
    end

    % Learning dates distance: 1 date, 2 distance, 3 std
    distancesAll    = [climateDistance stdDistance];
    distance        = distancesAll(~cellfun('isempty',distancesAll(:,1)),:);
    distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
    distancesBest   = distancesSort(1:nbImages,1);
    distSorted      = distancesSort(1:nbImages,2);
    stdSorted       = distancesSort(1:nbImages,3);
    sortedDates{qd} = currentQDate;
    sortedData{qd}  = cellfun(@double,distancesBest);
    sortedDist{qd}  = cellfun(@double,distSorted);
    sortedStd{qd}   = cellfun(@double,stdSorted);
end

sortedDatesAll = [sortedDates sortedData sortedDist sortedStd];
sortedDates    = sortedDatesAll;

end
