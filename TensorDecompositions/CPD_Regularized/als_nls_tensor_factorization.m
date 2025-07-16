function [A, B, C, D, output] = als_nls_tensor_factorization(T, M, k, rho, max_iters)
    % ALS with Nonlinear Least Squares for Tensor Factorization
    % Inputs:
    %   T - dxdxd Tensor
    %   M - dxd Matrix
    %   k - Rank of factorization
    %   rho - Regularization parameter
    %   max_iters - Number of ALS iterations
    % Outputs:
    %   A, B, C - d x k factor matrices
    %   D - k x k core matrix

    d = size(T, 1);

    % Initialize A, B, C, D with small random values
    A = rand(d, k); B = rand(d, k); C = rand(d, k); D = rand(k, k);
    
    % Define optimization options for fmincon
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');
    output.relerr = zeros(1, max_iters);

    for iter = 1:max_iters
        % Update A using nonlinear least squares
        T1 = unfold(T, 1);
        K_BC = khatri_rao(C, B);
        A = nls_solve(T1, K_BC, A, options);

        % Update B using nonlinear least squares
        T2 = unfold(T, 2);
        K_CA = khatri_rao(C, A);
        B = nls_solve(T2 + rho * M * B * D, K_CA + rho * D * B' * B, B, options);

        % Update C using nonlinear least squares
        T3 = unfold(T, 3);
        K_BA = khatri_rao(B, A);
        C = nls_solve(T3, K_BA, C, options);

        % Update D using nonlinear least squares
        D = nls_solve(M, B' * B, D, options);

        output.relerr(iter) = (1/2) * (norm(T - reshape(A * khatri_rao(C, B)', size(T)), "fro")^2 + (rho/2) * norm(M - B * D * B', "fro")^2);
    end
end

function X = nls_solve(Y, Z, X0, options)
    % Solve Nonlinear Least Squares with Constraints using fmincon
    d = size(Y, 1);
    k = size(X0, 2);
    
    % Define the objective function: ||Y - XZ||_F^2
    obj_fun = @(X) norm(Y - X * Z, 'fro')^2;

    % Constraints: X >= 0, sum(X(:)) = 1
    Aeq = ones(1, d * k); beq = 1; % Sum-to-one constraint
    lb = zeros(d, k); % Non-negativity constraint

    % Solve using fmincon
    X = fmincon(@(X) obj_fun(reshape(X, d, k)), X0(:), [], [], Aeq, beq, lb(:), [], [], options);
    X = reshape(X, d, k);
end

function T_unfold = unfold(T, mode)
    % Unfold tensor T along the specified mode
    d = size(T, 1);
    switch mode
        case 1
            T_unfold = reshape(permute(T, [1, 2, 3]), d, []);
        case 2
            T_unfold = reshape(permute(T, [2, 1, 3]), d, []);
        case 3
            T_unfold = reshape(permute(T, [3, 1, 2]), d, []);
        otherwise
            error('Invalid mode');
    end
end
