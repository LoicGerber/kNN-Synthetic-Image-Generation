function kge = computeKGE(simulated, observed)
    % Remove any invalid data points
    validMask = (simulated ~= -999 & observed ~= -999);
    simulated = simulated(validMask);
    observed = observed(validMask);

    % Compute the mean of observed and simulated values
    meanObserved = mean(observed);
    meanSimulated = mean(simulated);
    
    % Compute the standard deviation of observed and simulated values
    stdObserved = std(observed);
    stdSimulated = std(simulated);
    
    % Compute the correlation coefficient
    correlation = corr(simulated, observed);
    
    % Compute the bias ratio
    bias = meanSimulated / meanObserved;
    
    % Compute the variability ratio
    variability = (stdSimulated / meanSimulated) / (stdObserved / meanObserved);
    
    % Compute the KGE
    kge = 1 - sqrt((correlation - 1)^2 + (bias - 1)^2 + (variability - 1)^2);
end