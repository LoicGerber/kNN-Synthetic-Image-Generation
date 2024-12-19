function kge = computeKGE(simulated, observed)

    % Flatten the arrays to ensure they are treated as 1D vectors
    simulated = simulated(:);
    simulated = simulated(~isnan(simulated));
    observed = observed(:);
    observed = observed(~isnan(observed));

    % Compute the mean of observed and simulated values
    meanObserved = mean(observed,'omitnan');
    meanSimulated = mean(simulated,'omitnan');
    
    % Compute the standard deviation of observed and simulated values
    stdObserved = std(observed,'omitnan');
    stdSimulated = std(simulated,'omitnan');
    
    % Compute the correlation coefficient
    correlation = corr(simulated, observed);
    
    % Compute the bias ratio
    bias = meanSimulated / meanObserved;
    
    % Compute the variability ratio
    variability = (stdSimulated / meanSimulated) / (stdObserved / meanObserved);
    
    % Compute the KGE
    kge = 1 - sqrt((correlation - 1)^2 + (bias - 1)^2 + (variability - 1)^2);
end

