function combinedMetric = computeHellingerSPEM(x, y, spemFunc, hellingerFunc, spemHelW)
    spemDist       = spemFunc(x, y);
    hellingDist    = hellingerFunc(x, y);
    combinedMetric = (spemHelW(1) * (1 - spemDist)) + (spemHelW(2) * hellingDist);
%     combinedMetric = ((1 - spemDist) / 2) + (hellingDist / 2);
end