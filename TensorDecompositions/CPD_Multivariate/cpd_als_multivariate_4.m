function [A, B, C, output] = cpd_als_multivariate_4(P, k, options)
    [d, dn, ~] = size(P);
    n = dn/d;
    
    A = rand(dn, k);
    B = rand(dn, k);
    C = rand(dn, k);

    P1 = zeros(d, dn*dn);
    P2 = zeros(d, dn*dn);
    P3 = zeros(d, dn*dn);

    for i = 1:n
        i1 = 1 + d*(i-1);
        i2 = d*i;

        i3 = 1 + d*d*(i-1);
        i4 = d*d*i;

        P1(:, i3:i4) = mode_n_matricization(P(:, i1:i2, i1:i2), 1);
        P2(:, i3:i4) = mode_n_matricization(P(:, i1:i2, i1:i2), 2);
        P3(:, i3:i4) = mode_n_matricization(P(:, i1:i2, i1:i2), 3);
    end

    output.relerr = ones(3, options.maxIter);

    for iter = 1:options.maxIter
        A_old = A;
        B_old = B;
        C_old = C;

        A = P1 * khatri_rao(C, B) * pinv((C'*C).*(B'*B));
        A = Stochasticize(A);
    
        B = P2 * khatri_rao(C, A) * pinv((C'*C).*(A'*A));
        B = Stochasticize(B);
        
        C = P3 * khatri_rao(B, A) * pinv((B'*B).*(A'*A));
        C = Stochasticize(C);

        output.relerr(:, iter) = [A - A_old; B - B_old; C - C_old];

        if norm(A - A_old, 'fro') < options.tol && ...
           norm(B - B_old, 'fro') < options.tol && ...
           norm(C - C_old, 'fro') < options.tol
            fprintf('Converged at iteration %d\n', iter);
            break;
        end
    end

    A = Stochasticize(A);
    B = Stochasticize(B);
    C = Stochasticize(C);

    output.relerr = output.relerr(1:iter);
    output.iters = 1:iter;
end
