function [T, B, pi, output] = cpd_als_multivariate_3(P, nk, options)
    [d, ~, ~, n] = size(P);
    k = nk/n;

    A = cell(n, 1); B = cell(n, 1); C = cell(n, 1);
    P1 = cell(n, 1); P2 = cell(n, 1); P3 = cell(n, 1);

    epsilon = 1e-1;
    T = Stochasticize(eye(k) + epsilon * rand(k, k));
    pi = findStationaryDistribution(T')';

    for i = 1:n
        A{i} = Stochasticize(ones(d,k) + epsilon * rand(d, k));
        B{i} = Stochasticize(ones(d,k) + epsilon * rand(d, k));
        C{i} = Stochasticize(ones(d,k) + epsilon * rand(d, k));

        P1{i} = mode_n_matricization(P(:, :, :, i), 1);
        P2{i} = mode_n_matricization(P(:, :, :, i), 2);
        P3{i} = mode_n_matricization(P(:, :, :, i), 3);
    end

    output.relerr = ones(1, options.maxIter);


    for iter = 1:options.maxIter
        T_old = T;

        for i = 1:n
            Ai = P1{i} * khatri_rao(C{i}, B{i}) * pinv((C{i}'*C{i}).*(B{i}'*B{i}));
            A{i} = options.alpha * Ai + (1-options.alpha) * (B{i} * diag(pi) * T' * inv(diag(T * pi)));
            A{i} = Stochasticize(A{i});
            
            B{i} = P2{i} * khatri_rao(C{i}, A{i}) * pinv((C{i}'*C{i}).*(A{i}'*A{i}));
            B{i} = Stochasticize(B{i});

            Ci = P3{i} * khatri_rao(B{i}, A{i}) * pinv((B{i}'*B{i}).*(A{i}'*A{i}));
            C{i} = options.alpha * Ci + (1-options.alpha) * (B{i} * T);
            C{i} = Stochasticize(C{i});
        end


        Ts = cell(n,1);
        Ts{1} = Stochasticize(pinv(B{1}) * C{1});
        T = Ts{1}/n;

        for i = 2:n
            Ts{i} = Stochasticize(pinv(B{i}) * C{i});
            Pbest = PermutationFit(T, Ts{i}, true);

            Ts{i} = Pbest * Ts{i} * Pbest';

            T = T + Ts{i}/n;

            A{i} = A{i} * Pbest';
            B{i} = B{i} * Pbest';
            C{i} = C{i} * Pbest';

            output.col_perms{i} = Pbest';
        end

        T = options.rho * T_old + (1-options.rho) * T;
        pi = findStationaryDistribution(T')';
        output.relerr(iter) = norm(T-T_old, "fro");
        output.Ts = Ts;

        if(output.relerr(iter) < 1e-16)
            break;
        end

    end

    output.relerr = output.relerr(1:iter);
end