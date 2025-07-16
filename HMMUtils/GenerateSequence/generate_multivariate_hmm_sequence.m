function [observations, states] = generate_multivariate_hmm_sequence(T, O, pi, seq_length)
    % Generates a sequence of binned observations from a Hidden Markov Model.
    %
    % Parameters:
    % T : matrix (k x k)
    %     Transition matrix, where k is the number of states.
    % O : matrix (d x k x n)
    %     Observation matrix, where d is the number of possible observations.
    % pi : vector (1 x k)
    %     Initial state distribution.
    % seq_length : int
    %     Length of the sequence to generate.
    %
    % Returns:
    % obs_sequence : vector (n x seq_length)
    %     A sequence of observed indices (binned data).
    % states: a vector (1, seq_length)
    %     A sequence of hidden states.
    
    [d, k] = size(O{1});
    n = length(O);

    observations = zeros(n, seq_length);
    states = zeros(1, seq_length);

    cumsum_pi = cumsum(pi);
    cumsum_T = cumsum(T, 1);
    cumsum_O = cellfun(@(o) cumsum(o, 1), O, 'UniformOutput', false);

    % Sample initial state from pi
    state = find(mnrnd(1, pi));
    
    for t = 1:seq_length
        states(t) = state;

        % Sample observation from O(state, :)
        for i = 1:n
            observations(i, t) = find(rand < cumsum_O{i}(:, state), 1);
        end
        
        % Transition to the next state
        state = find(rand < cumsum_T(:, state), 1);
    end
end