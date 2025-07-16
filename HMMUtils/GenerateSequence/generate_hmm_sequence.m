function [observations, states] = generate_hmm_sequence(T, O, pi, seq_length)
    % Generates a sequence of binned observations from a Hidden Markov Model.
    %
    % Parameters:
    % T : matrix (k x k)
    %     Transition matrix, where k is the number of states.
    % O : matrix (d x k)
    %     Observation matrix, where d is the number of possible observations.
    % pi : vector (1 x k)
    %     Initial state distribution.
    % seq_length : int
    %     Length of the sequence to generate.
    %
    % Returns:
    % obs_sequence : vector (1 x seq_length)
    %     A sequence of observed indices (binned data).
    
    [d, k] = size(O);
    observations = zeros(1, seq_length);
    states = zeros(1, seq_length);

    cumsum_pi = cumsum(pi);
    cumsum_T = cumsum(T, 1);
    cumsum_O = cumsum(O, 1);

    % Sample initial state from pi
    state = find(mnrnd(1, pi));
    
    for t = 1:seq_length
        states(t) = state;

        % Sample observation from O(state, :)
        observations(t) = find(rand < cumsum_O(:, state), 1);
        
        % Transition to the next state
        state = find(rand < cumsum_T(:, state), 1);
    end
end