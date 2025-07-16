function [A, B, C, output] = cpd_als_3d_hmm_initialization(P, k, T_init, O_init, options)

    [d, ~, ~] = size(P);
    
    pi_init = findStationaryDistribution(T_init')';
    
    A = O_init * diag(pi_init) * T_init' * inv(diag(T_init*pi_init));
    B = O_init;
    C = O_init * T_init;


    P1 = mode_n_matricization(P, 1);
    P2 = mode_n_matricization(P, 2);
    P3 = mode_n_matricization(P, 3);

    output.relerr = ones(1, options.maxIter);

    for iter = 1:options.maxIter
        B = P2 * khatri_rao(C, A) * pinv((C'*C).*(A'*A));
        B = Stochasticize(B);
        
        C = P3 * khatri_rao(B, A) * pinv((B'*B).*(A'*A));
        C = Stochasticize(C);

        A = P1 * khatri_rao(C, B) * pinv((C'*C).*(B'*B));
        A = Stochasticize(A);

        output.relerr(iter) = norm(P - cpdgen({A, B, C}), "fro");
    end

    output.relerr = output.relerr(1:iter);
end