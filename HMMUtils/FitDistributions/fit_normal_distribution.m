function [mu, sigma] = fit_normal_distribution(P)
[d, k] = size(P);
y = 1:d;
mu = y * P;
sigma = sqrt(sum(P .* (y' * ones(1, k) - ones(d, 1) * mu) .^ 2, 1));
end