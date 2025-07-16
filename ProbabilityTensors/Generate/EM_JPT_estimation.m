function P3_n = EM_JPT_estimation(seq, d, rho, epsilon, nIter)

num_samples = length(seq);

% Initialize P3_n using empirical frequency counts
P3_n = zeros(d, d, d);

% Step 1: Compute raw counts of triplets from data
for t = 1:num_samples-2
    P3_n(seq(t), seq(t+1), seq(t+2)) = P3_n(seq(t), seq(t+1), seq(t+2)) + 1;
end

% Normalize to get initial probability estimates
P3_n = (P3_n + epsilon) / (sum(P3_n(:)) + numel(P3_n) * epsilon);

% EM Iteration
for iter = 1:nIter
    % E-Step: Compute expected counts by smoothing P3_n
    P3_smooth = (P3_n + epsilon) / (sum(P3_n(:)) + numel(P3_n) * epsilon);

    % M-Step: Update P3_n using the smoothed probability
    new_P3_n = zeros(size(P3_n));

    for t = 1:num_samples-2
        x = seq(t);
        y = seq(t+1);
        z = seq(t+2);

        % Compute expected contribution from each observed triplet
        new_P3_n(x, y, z) = rho * new_P3_n(x, y, z) + (1-rho) * P3_smooth(x, y, z);
    end

    % Normalize updated P3_n
    P3_n = (new_P3_n + epsilon) / (sum(new_P3_n(:)) + numel(new_P3_n) * epsilon);
    
end

end
