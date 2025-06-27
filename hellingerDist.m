function distanceOut = hellingerDist(x, y)
    if ndims(x) == 3 && ndims(y) == 3
        % Handle 3D arrays: compute Hellinger distance per layer
        numLayers = size(x, 3);
        distanceOut = inf(numLayers, 1);
        for n = 1:numLayers
            distanceOut(n) = hellingerDist2D(x(:,:,n), y(:,:,n));
        end
    else
        % Handle 2D input
        distanceOut = hellingerDist2D(x, y);
    end
end

function distance = hellingerDist2D(x, y)
    % Flatten and remove NaNs
    x = x(:);
    y = y(:);
    x = x(~isnan(x));
    y = y(~isnan(y));

    % Return maximum distance if inputs are empty or invalid
    minV = min([x; y]);
    maxV = max([x; y]);
    if isempty(x) || isempty(y) || minV >= maxV
        distance = 1;
        return;
    end

    % Define number of bins and edges
    numBins = 10;
    edges = linspace(minV, maxV, numBins + 1);

    % Compute and normalize histograms
    histX = histcounts(x, edges);
    histY = histcounts(y, edges);

    % Normalize to probability distributions
    histX = histX / sum(histX);
    histY = histY / sum(histY);

    % Compute Hellinger distance
    sqrtDiff = sqrt(histX) - sqrt(histY);
    distance = sqrt(sum(sqrtDiff.^2)) / sqrt(2);
end
