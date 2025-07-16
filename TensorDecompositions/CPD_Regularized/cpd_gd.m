function [A, B, C, D, output] = cpd_gd(T, M, k, rho, nu, max_iters)
    % ALS for Tensor Factorization with Matrix Regularization
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
    
    % Initialize A, B, C, D with random non-negative values
    A = rand(d, k); A = A ./ sum(A, 'all');
    B = rand(d, k); B = B ./ sum(B, 'all');
    C = rand(d, k); C = C ./ sum(C, 'all');
    D = rand(k, k); D = D ./ sum(D, 'all');

    T1 = unfold(T, 1);
    T2 = unfold(T, 2);
    T3 = unfold(T, 3);
    

    output.relerr = ones(1, max_iters);

    for iter = 1:max_iters
        % Update A
        dJdA = A*((C'*C).*(B'*B)) - T1 * khatri_rao(C, B);
        A = A - nu * dJdA;
        A = max(A, 0);
        A = A ./ sum(A, 'all');
        
        % Update B
        dJdB = B *((A'*A) .* (C'*C)) - T2*khatri_rao(C, A) + rho*(B*D'*(B')*B*D + B*D*(B')*B*D' - M*B*D' - M'*B*D);
        B = B - nu * dJdB;
        B = max(B, 0);
        B = B ./ sum(B, 'all');
        
        % Update C
        dJdC = C * ((B'*B) ./ (A' * A)) - T3 * khatri_rao(B, A);
        C = C - nu * dJdC;
        C = max(C, 0);
        C = C ./ sum(C, 'all');
        
        % Update D
        dJdD = B'*B*D*(B')*B - B'*M*B;
        D = D  - nu * dJdD;
        D = max(D, 0);
        D = D ./ sum(D, 'all');


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