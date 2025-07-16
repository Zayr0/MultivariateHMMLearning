function [P, errors] = estimate_joint_prob(sequence, num_bins, true_P)
    % Estimate a 3-way joint probability tensor from a sequence and track error relative to the true tensor.
    %
    % Inputs:
    %   sequence - A vector of binned data values (1 to num_bins)
    %   num_bins - The number of bins (maximum value in sequence)
    %   true_P - The true joint probability tensor for comparison
    %
    % Output:
    %   P - A num_bins x num_bins x num_bins tensor of estimated probabilities
    %   errors - A vector tracking error between estimated and true tensor

    % Initialize the probability tensor
    P = zeros(num_bins, num_bins, num_bins);

    s = 100;
    errors = ones(1, s);

    counter = 1;
    xerror = logspace(0, log10(length(sequence)-2), s);
    
    % Count occurrences of 3-wide windows
    for i = 1:length(sequence)-2
        a = sequence(i);
        b = sequence(i+1);
        c = sequence(i+2);
        P(a, b, c) = P(a, b, c) + 1;

        if(i >= xerror(counter))
            errors(counter) = norm(P / sum(P(:)) - true_P, 'fro');
            counter = counter +1;
        end
    end
    
    % Normalize to get probabilities
    P = P / sum(P(:));
end
