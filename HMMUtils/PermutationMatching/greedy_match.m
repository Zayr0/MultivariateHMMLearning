function [A_perm, row_perm, col_perm] = greedy_match(A, B)
    % Greedy matching for permutation of rows and columns
    % Inputs: 
    %   A - estimated k x k transition matrix
    %   B - true k x k transition matrix
    % Outputs:
    %   A_perm - permuted matrix A to match B
    %   row_perm - row permutation applied
    %   col_perm - column permutation applied

    k = size(A, 1); % Matrix size
    row_perm = zeros(1, k); 
    col_perm = zeros(1, k); 

    % ---- ROW MATCHING ----
    row_cost = zeros(k, k);
    for i = 1:k
        for j = 1:k
            row_cost(i, j) = norm(A(i, :) - B(j, :), 2); % Euclidean norm
        end
    end

    % Greedy assignment based on row distances
    available_rows = 1:k;
    available_targets = 1:k;
    for i = 1:k
        [min_val, idx] = min(row_cost(available_rows, available_targets), [], 'all', 'linear');
        [best_row, best_target] = ind2sub([numel(available_rows), numel(available_targets)], idx);
        
        row_perm(available_rows(best_row)) = available_targets(best_target);
        
        % Remove assigned row and target from future consideration
        available_rows(best_row) = [];
        available_targets(best_target) = [];
    end

    A = A(row_perm, :); % Apply row permutation

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
