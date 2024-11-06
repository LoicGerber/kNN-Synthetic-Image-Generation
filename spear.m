function rs = spear(X, Y)
    % Check that X and Y have the same length
    if length(X) ~= length(Y)
        error('X and Y must have the same length.');
    end
    
    % Compute the ranks of X and Y
    rank_X = tiedrank(X);
    rank_Y = tiedrank(Y);
    
    % Calculate the differences in ranks
    d = rank_X - rank_Y;
    
    % Compute Spearman's rank correlation coefficient
    n = length(X);
    rs = 1 - (6 * sum(d .^ 2)) / (n * (n^2 - 1));
end
