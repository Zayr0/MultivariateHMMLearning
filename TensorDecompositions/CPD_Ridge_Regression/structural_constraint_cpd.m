function [A, B, C, output] = structural_constraint_cpd(P3, k, rho, max_iters)
    % ALS for Tensor Factorization with structural constraints
    % Inputs:
    %   P3 - dxdxd Tensor
    %   k - Rank of factorization
    % Outputs:
    %   A, B, C - d x k factor matrices
    
    d = size(P3, 1);
    
    O = rand(d, k); O = Stochasticize(O);
    T = rand(k, k); T = Stochasticize(T);
    pi = findStationaryDistribution(T');

    A = Stochasticize(O * diag(pi) * T' / inv(diag(T * pi')));
    B = O;
    C = Stochasticize(O * T);

    T1 = unfold(P3, 1);
    T2 = unfold(P3, 2);
    T3 = unfold(P3, 3);

    output.relerr = ones(1, max_iters);

    for iter = 1:max_iters

        % Update A
        K_BC = khatri_rao(C, B);
        A = (T1 * K_BC + rho * (B * T')) * pinv(K_BC' * K_BC + rho*(T*T'));
        A = Stochasticize(abs(A));
        
        % Update B
        K_CA = khatri_rao(C, A);
        Binv = pinv(B);
        B = (T2 * K_CA + rho * (A*T + C*T' + Binv' * T)) * pinv(K_CA' * K_CA + rho*(eye(k) + T*T' + Binv * Binv'));
        B = Stochasticize(abs(B));
        
        % Update C
        K_BA = khatri_rao(B, A);
        Binv = pinv(B);
        C = (T3 * K_BA + rho*(B*T + Binv'*T)) * pinv(K_BA' * K_BA + rho * (eye(k) + Binv*Binv'));
        C = Stochasticize(abs(C));

        T = Stochasticize(Binv * C);
        
        output.relerr(iter) = (1/2) * norm(P3 - reshape(A * khatri_rao(C, B)', size(P3)), "fro")^2;
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