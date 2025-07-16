function [T, B, pi, output] = cpd_als_multivariate_2(P, nk, options)
    d = size(P, 1);
    n = size(P, 4);
    k = nk/n;

    A = cell(1, n); B = cell(1, n); C = cell(1, n);
    P1 = cell(1, n); P2 = cell(1, n); P3 = cell(1, n);

    Ts = cell(1, n);
    T = zeros(k,k);
    pis = cell(1, n);

    for i = 1:n
        A{i} = rand(d, k);
        B{i} = rand(d, k);
        C{i} = rand(d, k);

        P1{i} = mode_n_matricization(P(:, :, :, i), 1);
        P2{i} = mode_n_matricization(P(:, :, :, i), 2);
        P3{i} = mode_n_matricization(P(:, :, :, i), 3);
    end

    output.relerr = ones(1, options.maxIter);


    for iter = 1:options.maxIter
        T_old = T;
        

        % Fix T and solve all n cpd problems, obtain
        for i = 1:n
            A{i} = P1{i} * khatri_rao(C{i}, B{i}) * pinv((C{i}'*C{i}).*(B{i}'*B{i}));
            A{i} = Stochasticize(A{i});
        end

        for i = 1:n
            B{i} = P2{i} * khatri_rao(C{i}, A{i}) * pinv((C{i}'*C{i}).*(A{i}'*A{i}));
            B{i} = Stochasticize(B{i});
        end

        for i = 1:n
            C{i} = P3{i} * khatri_rao(B{i}, A{i}) * pinv((B{i}'*B{i}).*(A{i}'*A{i}));
            C{i} = Stochasticize(C{i});
        end
            
        T = zeros(k,k);
        for i = 1:n
            T = T + Stochasticize(pinv(B{i}) * C{i}) / n;
        end

        T = Stochasticize(T);
        pi = findStationaryDistribution(T');

        for i = 1:n
                A{i} = Stochasticize(B{i} * diag(pi) * T' / diag(T*pi'));
                C{i} = Stochasticize(B{i} * T);
        end

        output.relerr(iter) = norm(T-T_old, "fro");
    end

    output.relerr = output.relerr(1:iter);
end