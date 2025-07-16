function tt_cores = TT_ALS(tensor, ranks, max_iter, tol)
    % TT-ALS: Tensor Train Alternating Least Squares
    %
    % Inputs:
    % tensor  - The input tensor to approximate (N-dimensional array)
    % ranks   - TT-rank (vector of size N+1, where first and last are 1)
    % max_iter - Maximum number of ALS iterations (default: 50)
    % tol     - Convergence tolerance (default: 1e-6)
    %
    % Output:
    % tt_cores - Cell array containing TT-cores
    
    if nargin < 3, max_iter = 50; end
    if nargin < 4, tol = 1e-6; end
    
    dims = size(tensor);
    N = numel(dims);
    
    % Initialize TT-cores randomly
    tt_cores = cell(N, 1);
    for n = 1:N
        tt_cores{n} = randn(ranks(n), dims(n), ranks(n+1));
    end
    
    % Alternating Least Squares
    for iter = 1:max_iter
        prev_cores = tt_cores;
        
        for n = 1:N
            % Compute unfolding matrix of the tensor (fixing core n)
            unfolding = reshape(tensor, [prod(dims(1:n)), prod(dims(n+1:end))]);
            
            % Compute the left and right environment matrices
            left_env = eye(ranks(n));
            right_env = eye(ranks(n+1));
            
            if n > 1
                left_env = contract_left(tt_cores, n);
            end
            if n < N
                right_env = contract_right(tt_cores, n);
            end
            
            % Solve least squares for core update
            core_update = left_env \ unfolding / right_env';
            tt_cores{n} = reshape(core_update, [ranks(n), dims(n), ranks(n+1)]);
        end
        
        % Check for convergence
        diff = max(cellfun(@(c1, c2) norm(c1(:) - c2(:)), tt_cores, prev_cores));
        if diff < tol
            break;
        end
    end
end

function left_env = contract_left(cores, n)
    % Computes the left environment for core n
    left_env = eye(size(cores{1}, 1));
    for i = 1:n-1
        left_env = left_env * unfold(cores{i}, 3);
    end
end

function right_env = contract_right(cores, n)
    % Computes the right environment for core n
    right_env = eye(size(cores{end}, 3));
    for i = numel(cores):-1:n+1
        right_env = unfold(cores{i}, 1) * right_env;
    end
end

function mat = unfold(core, mode)
    % Unfolds a TT-core along the specified mode
    dims = size(core);
    if mode == 1
        mat = reshape(core, dims(1), []);
    elseif mode == 3
        mat = reshape(core, [], dims(3));
    else
        mat = reshape(core, [dims(1) * dims(2), dims(3)]);
    end
end
