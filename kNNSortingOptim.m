function sortedDate = kNNSortingOptim(sortedDates,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW,Pre_ShortW,Tavg_LongW,Tmin_LongW,Tmax_LongW,Pre_LongW,spemW,helW,metricKNN,nbImages,saveOptimPrep,inputDir)

%
%
%
% REDO DOCUMENTATION
%
%
%

% 1 Q date, 2 L dates, 3 climate distances
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
weightsShort  = [Pre_ShortW,Tavg_ShortW,Tmin_ShortW,Tmax_ShortW];
weightsLong   = [Pre_LongW,Tavg_LongW,Tmin_LongW,Tmax_LongW];
% if ~isempty(addVars)
%     error('addVars weights handling not implemented yet')
%     %idxAddVars     = contains(weightsNames,addVars);
%     %weightsAddVars = bayesWeights(:,idxAddVars);
% else
%     weightsAddVars = [];
% end

totQDates   = size(distances,1);

sortedDate  = cell(totQDates, 1);
sortedData  = cell(totQDates, 1);
sortedDist  = cell(totQDates, 1);

for qd = 1:totQDates
    climateDistance = cell(totQDates, 2);
    currentQDate = distances{qd,1};
    totLDates    = size(distances{qd,2},1);
    for ld = 1:totLDates
        currentLDate = distances{qd,2}{ld,1};
        
        % Climate distance
        climateDistance{ld,1} = currentLDate;
        climateDistance{ld,2} = distances{qd,3}{ld,1};
        if metricKNN == 5
            climateDistance{ld,2}(1,:) = distances{qd,3}{ld,1}(1,:) .* weightsShort;
            climateDistance{ld,2}(2,:) = distances{qd,3}{ld,1}(2,:) .* weightsLong;
            climateDistance{ld,3}(1,:) = distances{qd,3}{ld,2}(1,:) .* weightsShort;
            climateDistance{ld,3}(2,:) = distances{qd,3}{ld,2}(2,:) .* weightsLong;

            climateDistance{ld,2}(1,:) = (climateDistance{ld,2}(1,:) * spemW) + (climateDistance{ld,3}(1,:) * helW);
            climateDistance{ld,2}(2,:) = (climateDistance{ld,2}(2,:) * spemW) + (climateDistance{ld,3}(2,:) * helW);
        else
            climateDistance{ld,2}(1,:) = climateDistance{ld,2}(1,:) .* weightsShort;
            climateDistance{ld,2}(2,:) = climateDistance{ld,2}(2,:) .* weightsLong;
        end
        climateDistance{ld,2} = sum(sum(climateDistance{ld,2},1),2);
    end

    % Learning dates distance: 1 date, 2 distance, 3 std
    distance       = climateDistance(~cellfun('isempty',climateDistance(:,1)),:);
    distancesSort  = sortrows(distance,2); % Sort rows in ascending order according to column 2
    distancesBest  = distancesSort(1:nbImages,1);
    distSorted     = distancesSort(1:nbImages,2);
    sortedDate{qd} = currentQDate;
    sortedData{qd} = cell2mat(distancesBest);
    sortedDist{qd} = cell2mat(distSorted);
end

sortedDatesAll = [sortedDate sortedData sortedDist];
sortedDate     = sortedDatesAll;

end
