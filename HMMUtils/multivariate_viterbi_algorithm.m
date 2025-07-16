function [vPath, delta] = multivariate_viterbi_algorithm(obs, T, Os, pi, alpha, mvWindowSize)
    % obs: n × m matrix (n sequences, m time steps) Note that these values
    % are discretized to be between 1 and d = # of symbols/emissions
    % T: k × k transition matrix
    % Os: cell array of n emission matrices (d × k each), where d = # of symbols
    % pi: k × 1 initial state probabilities
    % alpha: scalar in [0,1] controlling emission vs transition emphasis

    if nargin < 5
        alpha = 1; % default: full emission-based Viterbi
    elseif nargin < 6
        mvWindowSize = 1;
    end

    logeps = -1e10;

    [n, m] = size(obs);               % n: # sequences, m: sequence length
    [d, k] = size(Os{1});             % d: emission symbols per variable, k: states

    % Validate inputs
    if (length(Os) ~= n || size(T, 1) ~= k || size(T, 2) ~= k || length(pi) ~= k)
        error("Dimension mismatch in inputs");
    end

    % Convert emission matrices from cell to array: Os(i, symbol, state)
    Omat = zeros(n, d, k);
    for i = 1:n
        Omat(i, :, :) = Os{i};
    end

    % Replace 0s to avoid log(0)
    T(T <= 0) = 1e-50;
    Omat(Omat <= 0) = 1e-50;
    pi(pi <= 0) = 1e-50;

    logT = log(T);
    logO = log(Omat);
    logPi = log(pi);

    delta = zeros(m, k);    % log-delta
    prev = zeros(m, k);     % traceback

    % --- Initialization ---
    for s = 1:k
        log_emission = 0;
        for i = 1:n
            symbol = obs(i, 1);
            log_emission = log_emission + logO(i, symbol, s);
        end
        delta(1, s) = 2 * (alpha * log_emission + (1 - alpha) * logPi(s));
    end
    delta(1, :) = log(exp(delta(1, :))./sum(exp(delta(1, :))));

    % --- Loop ---
    for t = 2:m         % Check each time step
        for s = 1:k     % Check each current state
            max_score = logeps;
            argmax_r = 0;

            % Emission log-probability for state s at time t
            log_emission = 0;
            for i = 1:n
                symbol = obs(i, t);
                log_emission = log_emission + logO(i, symbol, s);
            end

            for r = 1:k         % Check each previous state with respect to the current state
                log_transition = logT(s, r);
                score = 2 * (0.5 * delta(t-1, r) + (alpha * log_emission + (1 - alpha) * log_transition));

                if score > max_score
                    max_score = score;
                    argmax_r = r;
                end
            end

            delta(t, s) = max_score;
            % delta(t, :) = log(exp(delta(t, :))./sum(exp(delta(t, :))));

            prev(t, s) = argmax_r;
        end
        % maxLog = max(delta(t, :));  % for numerical stability
        % logSumExp = maxLog + log(sum(exp(delta(t, :) - maxLog)));
        % delta(t, :) = delta(t, :) - logSumExp;
    end

    % ensure prev has all the same entries
    % prev = mode(prev, 2) * ones(1, k);

    % --- Backtracking ---
    vPath = zeros(m, 1);
    [~, vPath(m)] = max(delta(m, :));
    for t = m-1:-1:1
        vPath(t) = prev(t+1, vPath(t+1));
    end
    
    if(mvWindowSize > 0)
        vPath = movmode(vPath, mvWindowSize);
    end
end


function smoothed = movmode(path, window_size)
    % Ensure window size is odd
    if mod(window_size, 2) == 0
        window_size = window_size + 1;
    end

    half_w = floor(window_size / 2);
    n = length(path);
    smoothed = zeros(size(path));

    for i = 1:n
        % Compute window indices with edge clamping
        start_idx = max(1, i - half_w);
        end_idx = min(n, i + half_w);

        window = path(start_idx:end_idx);
        smoothed(i) = mode(window);
    end
end

% function corrected_matrix = fill_majority_rows(A)
%     [m, k] = size(A);
%     corrected_matrix = zeros(m, k);  % Preallocate for speed
% 
%     for i = 1:m
%         row = A(i, :);
%         % Get unique values and their counts
%         [unique_vals, ~, idx] = unique(row);
%         counts = accumarray(idx', 1);
%         % Find value with the maximum count (mode)
%         [~, max_idx] = max(counts, [],'all');
%         majority_val = unique_vals(max_idx);
%         % Fill the row with the majority value
%         corrected_matrix(i, :) = majority_val * ones(1, k);
%     end
% end



%% Multivariate Viterbi without LOG - Might have numerical errors

% function [vPath, delta] = multivariate_viterbi_algorithm(obs, T, Os, pi)
%     transfrac = 1;
% 
%     [n, m] = size(obs);
%     [d, k] = size(Os{1});
% 
%     if (length(Os) ~= n || size(T, 1) ~= k || size(T, 2) ~= k || size(pi, 1) ~= k)
%         warning("Some dimensions are wrong");
%     end
% 
%     if(iscell(Os))
%         Omat = zeros(n, d, k);
%         for i = 1:n
%             Omat(i, :, :) = Os{i};
%         end
%         Os = Omat;
%     end
%     T(T <= 0) = 1e-50;
%     Os(Os <= 0) = 1e-50;
%     pi(pi <= 0) = 1e-50;
% 
%     delta = zeros(m, k);
%     prev = zeros(m, k);
% 
%     for s = 1:k
%         delta(1, s) = pi(s);
% 
%         for i = 1:n
%             symbol = obs(1, i);
%             delta(1, s) = delta(1, s) * Os(i, symbol, s);
%         end
%     end
%     delta(1, :) = delta(1, :) ./ sum(delta(1, :));
% 
% 
%     for t = 2:m-1             % Check each timestep
%         for s = 1:k         %Check each current state
%             max_prob = 0;
%             argmax_prev = 0;
%             for r = 1:k     %Check each previous state
%                  trans_prob = transfrac * delta(t-1, r) * T(s, r);
%                  emission_prob = 1;
%                  for i = 1:n
%                      symbol = obs(i, t);
%                      emission_prob = emission_prob * Os(i, symbol, s);
%                  end
%                  new_prob = trans_prob * emission_prob;
%                  if(new_prob > max_prob)
%                     max_prob = new_prob;
%                     argmax_prev = r;
%                  end
%             end
%             delta(t, s) = max_prob;
%             prev(t, s) = argmax_prev;
%         end
%         delta(t, :) = delta(t, :) ./ sum(delta(t, :));
%     end
% 
%     [~, vPath] = max(delta, [], 2);
% end