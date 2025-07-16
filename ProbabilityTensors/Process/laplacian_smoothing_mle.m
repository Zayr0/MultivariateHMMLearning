function P3_mle = laplacian_smoothing_mle(P3_n, alpha)
    % ESTIMATE_JOINT_PROB_MLE - Corrects an empirical joint probability tensor
    % using Maximum Likelihood Estimation (MLE) with Laplace smoothing.
    %
    % Inputs:
    %   P3_n  - Empirical joint probability tensor (d x d x d) (estimated from samples)
    %   alpha - Smoothing parameter (default: 1 for Laplace smoothing)
    %
    % Output:
    %   P3_mle - Adjusted joint probability tensor (d x d x d)
    
    % Set default smoothing parameter if not provided
    if nargin < 2
        alpha = 1; % Laplace smoothing
    end
    
    % Get dimensions of the tensor
    d = size(P3_n, 1);
    
    % Compute the total number of observed triplets (assuming P3_n is unnormalized)
    total_count = sum(P3_n(:));
    
    % Compute the MLE estimate with Laplace smoothing
    P3_mle = (P3_n + alpha) / (total_count + alpha * d^3);
    
    % Ensure numerical stability (avoid rounding errors)
    P3_mle = max(P3_mle, 1e-10); % Prevent zero probabilities
    P3_mle = P3_mle / sum(P3_mle(:)); % Renormalize
    
end