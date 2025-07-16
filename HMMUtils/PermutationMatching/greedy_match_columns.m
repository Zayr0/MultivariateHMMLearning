function [A_perm, col_perm] = greedy_match_columns(A, B)
    % Greedy matching for permutation of rows and columns
    % Inputs: 
    %   A - estimated k x k transition matrix
    %   B - true k x k transition matrix
    % Outputs:
    %   A_perm - permuted matrix A to match B
    %   row_perm - row permutation applied
    %   col_perm - column permutation applied

    [d, k] = size(A);
    col_perm = zeros(1, k); 


    % ---- COLUMN MATCHING ----
    col_cost = zeros(k, k);
    for i = 1:k
        for j = 1:k
            col_cost(i, j) = norm(A(:, i) - B(:, j), 2); % Euclidean norm
        end
    end

    % Greedy assignment based on column distances
    available_cols = 1:k;
    available_targets = 1:k;
    for i = 1:k
        [min_val, idx] = min(col_cost(available_cols, available_targets), [], 'all', 'linear');
        [best_col, best_target] = ind2sub([numel(available_cols), numel(available_targets)], idx);
        
        col_perm(available_cols(best_col)) = available_targets(best_target);
        
        % Remove assigned column and target from future consideration
        available_cols(best_col) = [];
        available_targets(best_target) = [];
    end

    A_perm = A(:, col_perm); % Apply column permutation
end
