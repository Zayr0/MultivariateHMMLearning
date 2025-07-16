function contrast = matrix_contrast(A)
    mu = mean(A(:));  % Compute the mean of all elements
    sigma = std(A(:)); % Compute the standard deviation
    contrast = sigma / mu; % Coefficient of variation
end