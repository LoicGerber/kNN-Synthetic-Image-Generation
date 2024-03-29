function sortedDates = kNNSortingOptim(sortedDates,addVars,Et_W,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,nbImages,saveOptimPrep,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% 1 Q date, 2 L dates, 3 target distance, 4 addVars distances, 5 climate distances, 6 std
if saveOptimPrep == true
    distances = load(fullfile(inputDir,'KNNDistances.mat'));
    distances = distances.sortedDates;
else
    distances = sortedDates;
end

% Assign different weights
% weightsNames  = {bayesWeights.Name};
% idxTarget     = contains(weightsNames,strcat(var,'_W'));
% weightsTarget = bayesWeights(:,idxTarget);
% idxShort      = contains(weightsNames,'Short');
% weightsShort  = bayesWeights(:,idxShort);
% idxLong       = contains(weightsNames,'Long');
% weightsLong   = bayesWeights(:,idxLong);
weightsTarget = Et_W;
weightsShort  = [Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW];
weightsLong   = [Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW];
if ~isempty(addVars)
    error('addVars weights handling not implemented yet')
    %idxAddVars     = contains(weightsNames,addVars);
    %weightsAddVars = bayesWeights(:,idxAddVars);
else
    weightsAddVars = [];
end

totQDates   = size(distances,1);
totLDates   = size(distances{1,2},1);

targetDistance  = cell(totQDates, 1);
addVarsDistance = cell(totQDates, 2);
climateDistance = cell(totQDates, 2);

sortedDates = cell(totQDates, 1);
sortedData  = cell(totQDates, 1);
sortedDist  = cell(totQDates, 1);

for qd = 1:totQDates
    currentQDate = distances{qd,1};
    for ld = 1:totLDates
        currentLDate = distances{1,2}{ld,1};
        % Target variable comparison
        targetDistance{ld,1} = distances{qd,3}{ld,1};
        targetDistance{ld,1} = targetDistance{ld,1} .* weightsTarget;

        % Additional variable comparison
        if ~isempty(addVars)
            if numel(addVars) == 1
                addVarsDistance{ld,1} = distances{qd,4}{ld,1};
                addVarsDistance{ld,1} = addVarsDistance{ld,1} .* weightsAddVars;
                addVarsDistance{ld,2} = addVarsDistance{ld,2} .* weightsAddVars;
            else
                addVarsDistance{:,ld} = distances{:,ld};
                addVarsDistance{:,ld} = addVarsDistance{:,ld} .* weightsAddVars;
                addVarsDistance{ld,2} = addVarsDistance{ld,2} .* weightsAddVars;
            end
        else
            addVarsDistance{ld,1} = 0;
            addVarsDistance{ld,2} = 0;
        end

        % Climate distance
        climateDistance{ld,1} = currentLDate;
        climateDistance{ld,2} = distances{qd,5}{ld,1};
        climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* weightsShort;
        climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* weightsLong;
        climateDistance{ld,2} = sum(climateDistance{ld,2},1,'omitnan');
        climateDistance{ld,2} = sum(climateDistance{ld,2},2,'omitnan')+targetDistance{ld,1}+addVarsDistance{ld,1};
    end

    % Learning dates distance: 1 date, 2 distance, 3 std
    distance        = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
    distancesSort   = sortrows(distance,2); % Sort rows in ascending order according to column 2
    distancesBest   = distancesSort(1:nbImages,1);
    distSorted      = distancesSort(1:nbImages,2);
    sortedDates{qd} = currentQDate;
    sortedData{qd}  = cell2mat(distancesBest);
    sortedDist{qd}  = cell2mat(distSorted);
end

sortedDatesAll = [sortedDates sortedData sortedDist];
sortedDates    = sortedDatesAll;

end
