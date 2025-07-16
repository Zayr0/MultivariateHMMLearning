function [mu, Sigma, triplets] = fit_joint_gaussian_triplets(data)
% Fits a 3-way joint Gaussian distribution to sliding triplets in a 1D data vector
% 
% Inputs:
%   data    - vector of length m
%
% Outputs:
%   mu      - 1x3 vector of means
%   Sigma   - 3x3 covariance matrix
%   triplets - (m-2)x3 matrix of consecutive triplets

    % Input validation
    if ~isvector(data)
        error('Input must be a 1D vector.');
    end
    if length(data) < 3
        error('Input must contain at least 3 elements.');
    end

    data = data(:);  % ensure column vector
    m = length(data);
    
    % Form overlapping triplets: [x1 x2 x3; x2 x3 x4; ...]
    triplets = [data(1:end-2), data(2:end-1), data(3:end)];
    
    % Fit multivariate normal distribution
    mu = mean(triplets);         % 1x3 mean vector
    Sigma = cov(triplets);       % 3x3 covariance matrix
end
