function [alpha, beta] = fit_beta_distribution(P)
    [d, k] = size(P);
    y = linspace(0, 1, d);
    mu = y * P;
    sigma = sqrt(sum(P .* (y' * ones(1, k) - ones(d, 1) * mu) .^ 2, 1));
    
    alpha = abs(mu .* ((mu .* (ones(1, k) - mu) ./ sigma) - ones(1, k)));
    beta = abs((ones(1, k) - mu) .* ((mu .* (ones(1, k) - mu) ./ sigma) - ones(1, k)));
end
