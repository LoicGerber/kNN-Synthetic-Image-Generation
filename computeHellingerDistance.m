% Function to calculate Hellinger distance based on histograms
function distance = computeHellingerDistance(values1, values2)
    % Create histograms with automatic number of bins
    %numBins = max(size(histcounts(values1),2), size(histcounts(values2),2));
    % Create histograms with fixed number of bins
    if min(min([values1; values2])) < max(max([values1; values2]))
        numBins = 10;
        edges = linspace(min(min([values1; values2])), max(max([values1; values2])), numBins + 1);
    
        % Compute histograms
        hist1 = histcounts(values1, edges);
        hist2 = histcounts(values2, edges);
    
        % Normalize histograms
        normHist1 = hist1 / sum(hist1);
        normHist2 = hist2 / sum(hist2);
    
        % Compute Hellinger distance
        distance = 1 / sqrt(2) * sqrt(sum((sqrt(normHist1) - sqrt(normHist2)).^2));
    else
        distance = 0;
    end
end