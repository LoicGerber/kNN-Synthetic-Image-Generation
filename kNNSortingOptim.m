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
    distances = load(fullfile(inputDir, 'KNNDistances.mat'));
    distances = distances.sortedDates;
else
    distances = sortedDates;
end

weightsShort = [Pre_ShortW, Tavg_ShortW, Tmin_ShortW, Tmax_ShortW];
weightsLong  = [Pre_LongW, Tavg_LongW, Tmin_LongW, Tmax_LongW];

totQDates = size(distances, 1);

sortedDate = cell(totQDates, 1);
sortedData = cell(totQDates, 1);
sortedDist = cell(totQDates, 1);

for qd = 1:totQDates
    currentQDate = distances{qd, 1};
    learningDates = cell2mat(distances{qd, 2});
    distCell = distances{qd, 3};

    totLDates = numel(learningDates);
    climateDistance = NaN(totLDates, 2);

    % Fill first column with learning dates
    climateDistance(:, 1) = learningDates;

    % Preallocate and concatenate distance matrices
    cdist = cat(3, distCell{:, 1});
    hdist = cat(3, distCell{:, 2});

    % Apply weights and aggregate
    D1 = squeeze(cdist(1, :, :))';
    D2 = squeeze(cdist(2, :, :))';

    D1_weighted = D1 .* weightsShort;
    D2_weighted = D2 .* weightsLong;

    if metricKNN == 5
        H1 = squeeze(hdist(1, :, :))';
        H2 = squeeze(hdist(2, :, :))';

        D1_weighted = D1_weighted * spemW + (H1 .* weightsShort) * helW;
        D2_weighted = D2_weighted * spemW + (H2 .* weightsLong)  * helW;
    end

    totalDistance = sum(D1_weighted + D2_weighted, 2);

    climateDistance(:, 2) = totalDistance;

    % Sort and select
    distancesSort = sortrows(climateDistance, 2);
    distancesBest = distancesSort(1:nbImages, 1);
    distSorted    = distancesSort(1:nbImages, 2);

    sortedDate{qd} = currentQDate;
    sortedData{qd} = distancesBest;
    sortedDist{qd} = distSorted;
end

sortedDatesAll = [sortedDate, sortedData, sortedDist];
sortedDate = sortedDatesAll;

end
