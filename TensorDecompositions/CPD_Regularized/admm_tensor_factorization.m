function [A, B, C, D, output] = admm_tensor_factorization(T, M, k, rho, gamma, max_iters)
    % ADMM for Tensor Factorization with Matrix Regularization
    % Inputs:
    %   T - dxdxd Tensor
    %   M - dxd Matrix
    %   k - Rank of factorization
    %   rho - Regularization parameter
    %   gamma - ADMM penalty parameter
    %   max_iters - Number of ADMM iterations
    % Outputs:
    %   A, B, C - d x k factor matrices
    %   D - k x k core matrix

    d = size(T, 1);

    % Initialize variables
    A = rand(d, k); A = A ./ sum(A, 'all');
    B = rand(d, k); B = B ./ sum(B, 'all');
    C = rand(d, k); C = C ./ sum(C, 'all');
    D = rand(k, k); D = D ./ sum(D, 'all');

    Z_A = A; Z_B = B; Z_C = C; Z_D = D;
    lambda_A = zeros(d, k); lambda_B = zeros(d, k); lambda_C = zeros(d, k); lambda_D = zeros(k, k);

    output.relerr = zeros(1, max_iters);

    for iter = 1:max_iters
        % Update A
        T1 = unfold(T, 1);
        K_BC = khatri_rao(C, B);
        A = (T1 * K_BC + gamma * (Z_A - lambda_A)) * pinv(K_BC' * K_BC + gamma * eye(k));
        A = Stochasticize(abs(A));

        % Update B
        T2 = unfold(T, 2);
        K_CA = khatri_rao(C, A);
        B = (T2 * K_CA + rho * M * B * D + gamma * (Z_B - lambda_B)) * pinv(K_CA' * K_CA + rho * D * B' * B + gamma * eye(k));
        B = Stochasticize(abs(B));

        % Update C
        T3 = unfold(T, 3);
        K_BA = khatri_rao(B, A);
        C = (T3 * K_BA + gamma * (Z_C - lambda_C)) * pinv(K_BA' * K_BA + gamma * eye(k));
        C = Stochasticize(abs(C));

        % Update D
        D = pinv(B' * B + gamma * eye(k)) * (B' * M * B + gamma * (Z_D - lambda_D)) * pinv(B' * B + gamma * eye(k));
        D = Stochasticize(abs(D));

        % Update auxiliary variables (Projection step)
        Z_A = max(A + lambda_A, 0); Z_A = Z_A / sum(Z_A, 'all');
        Z_B = max(B + lambda_B, 0); Z_B = Z_B / sum(Z_B, 'all');
        Z_C = max(C + lambda_C, 0); Z_C = Z_C / sum(Z_C, 'all');
        Z_D = max(D + lambda_D, 0); Z_D = Z_D / sum(Z_D, 'all');

        % Update Lagrange multipliers
        lambda_A = lambda_A + A - Z_A;
        lambda_B = lambda_B + B - Z_B;
        lambda_C = lambda_C + C - Z_C;
        lambda_D = lambda_D + D - Z_D;

        output.relerr(iter) = (1/2) * (norm(T - reshape(A * khatri_rao(C, B)', size(T)), "fro")^2 + (rho/2) * norm(M - B * D * B', "fro")^2);
    end
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
