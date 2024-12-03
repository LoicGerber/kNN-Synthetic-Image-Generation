% Function to calculate Hellinger distance based on histograms
function distance = computeHellingerDistance(x, y)
    % Remove non-positive values
%     values1 = values1(values1 > 0);
%     values2 = values2(values2 > 0);

    % Return maximum distance if inputs are empty or invalid
    minV = min([x(:); y(:)]);
    maxV = max([x(:); y(:)]);
    if isempty(x) || isempty(y) || minV >= maxV
        distance = 1;
        return;
    end

    % Define number of bins and edges
    numBins = 10;
    edges = linspace(minV, maxV, numBins + 1);

    % Compute and normalize histograms
    histX = histcounts(x, edges, 'Normalization', 'probability');
    histY = histcounts(y, edges, 'Normalization', 'probability');

    % Compute Hellinger distance
    sqrtDiff = sqrt(histX) - sqrt(histY);
    distance = sqrt(sum(sqrtDiff.^2)) / sqrt(2);
end

