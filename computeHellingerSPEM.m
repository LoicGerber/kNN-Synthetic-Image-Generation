function combinedMetric = computeHellingerSPEM(x, y, spemFunc, hellingerFunc)
    spemDist = spemFunc(x, y);
    hellingDist = hellingerFunc(x, y);
    combinedMetric = ((1 - spemDist) / 2) + (hellingDist / 2);
end