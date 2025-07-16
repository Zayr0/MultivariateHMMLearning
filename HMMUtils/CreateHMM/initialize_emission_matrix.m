function emission_matrix = initialize_emission_matrix(d, k, n, mean_range, std_range)
    % Initialize a d x k x n emission matrix where each state's emission probabilities follow a Gaussian distribution.
    % The Gaussians are spaced apart to minimize linear dependence.
    %
    % Parameters:
    % - d (int): Number of emissions
    % - k (int): Number of states
    % - n (int): Number of temporal variables (or separate emission matrices)
    % - mean_range (vector): Range for mean selection
    % - std_range (vector): Range for standard deviation selection
    %
    % Returns:
    % - emission_matrix: A (d, k, n) emission probability matrix

    alpha = 1e-9; % Small value to ensure numerical stability
    
    if nargin < 4
        mean_range = [0.1, 0.9] * d; % Ensure means are spread across the range
    end
    if nargin < 5
        std_range = [0.05, 0.15] * d; % Set standard deviation range
    end
    
    emission_matrix = cell(n, 1);
    
    for var = 1:n
        O = zeros(d, k);
        
        for state = 1:k
            % Randomly perturb the mean slightly
            mean_val = mean_range(1) + (mean_range(2) - mean_range(1)) * rand();
            
            % Randomly choose a standard deviation
            std_val = std_range(1) + (std_range(2) - std_range(1)) * rand();
            
            % Generate Gaussian distribution
            emissions = exp(-0.5 * ((1:d) - mean_val).^2 / std_val^2);
            emissions = emissions / sum(emissions); % Normalize
            
            % Store in matrix
            O(:, state) = emissions + alpha;
        end
        
        % Normalize the matrix to ensure it remains stochastic
        O = Stochasticize(O);
        emission_matrix{var} = O(:, randperm(k));

    end
    
    if n == 1
        emission_matrix = emission_matrix{1};
    elseif n < 1
        warning("wrong choice of n");
    end
end

