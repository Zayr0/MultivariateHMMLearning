function P_best = PermutationFitEmission(O, O_hat)
    k = size(O, 2);  % Get the number of states
    
    % Generate all possible permutations of indices [1:k]
    permutations = perms(1:k);
    num_perms = size(permutations, 1);
    
    best_error = inf;  % Initialize best error as a large value
    best_perm = 1:k;   % Default best permutation

    % Try all possible permutations
    for i = 1:num_perms
        perm = permutations(i, :);  % Get the current permutation
        P = eye(k);                 % Start with identity matrix
        P = P(perm, :);              % Permute rows
        
        % Apply permutation to O_hat
        Q_permuted = O_hat * P';  % Permute both rows and columns

        % Compute Frobenius norm ||O - Q_permuted||
        err = norm(O - Q_permuted, 'fro');

        % Store the best permutation
        if err < best_error
            best_error = err;
            best_perm = perm;
        end
    end

    % Apply the best permutation found
    P_best = eye(k);
    P_best = P_best(best_perm, :);
    
end