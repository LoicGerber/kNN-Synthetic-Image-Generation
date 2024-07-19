function nse = computeNSE(simulated, observed)
    % Remove any invalid data points
    validMask = (simulated ~= -999 & observed ~= -999);
    simulated = simulated(validMask);
    observed = observed(validMask);
    
    % Compute the mean of observed values
    meanObserved = mean(observed);
    
    % Compute the numerator (sum of squared differences between observed and simulated values)
    numerator = sum((observed - simulated).^2);
    
    % Compute the denominator (sum of squared differences between observed values and mean observed value)
    denominator = sum((observed - meanObserved).^2);
    
    % Compute the NSE
    nse = 1 - (numerator / denominator);
end
