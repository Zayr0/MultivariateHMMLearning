function [T, B, pi, output] = multivariate_coupled_cpd_3(P, K, D, N, alpha, rho, options)
    
    % Initialize the factor matrices and transition matrix
    A = cell(N, 1); B = cell(N, 1); C = cell(N, 1); B_prev = cell(N, 1);
    P1 = cell(N, 1); P2 = cell(N, 1); P3 = cell(N, 1);
    
    if isfield(options, 'epsilon')
        epsilon = options.epsilon;
    else
        epsilon = 1e-3;
    end
    
    T = Stochasticize(eye(K) + epsilon * rand(K, K));
    pi = findStationaryDistribution(T')';

    for n = 1:N
        A{n} = Stochasticize(ones(D,K) + epsilon * rand(D, K));
        B{n} = Stochasticize(ones(D,K) + epsilon * rand(D, K));
        C{n} = Stochasticize(ones(D,K) + epsilon * rand(D, K));

        P1{n} = mode_n_matricization(P{n}, 1);
        P2{n} = mode_n_matricization(P{n}, 2);
        P3{n} = mode_n_matricization(P{n}, 3);
    end
    B_prev = B;
    
    % Allocate memory for error metrics: convergence of transition matrix,
    % and tensor reconstruction error.
    output.tconv = ones(1, options.maxIter);
    output.recerr = ones(N, options.maxIter);
    output.rellrecerr = ones(N, options.maxIter);
    output.emiscorr = ones(N, options.maxIter);


    alpha = 1/(1 + alpha);


    for iter = 1:options.maxIter
        T_prev = T; % Save previous T to check convergence later
        
        % Calculate factor matrices of every joint probability tensor using
        % ALS with partial dependece on common T
        for n = 1:N
            An = P1{n} * khatri_rao(C{n}, B{n}) * pinv((C{n}'*C{n}).*(B{n}'*B{n}));
            A{n} = Stochasticize(alpha* An + (1-alpha) * (B{n} * diag(pi) * T' / inv(diag(T * pi))));
            
            Bn = P2{n} * khatri_rao(C{n}, A{n}) * pinv((C{n}'*C{n}).*(A{n}'*A{n}));
            B{n} = Stochasticize(alpha * Bn + (1-alpha) * B{n});

            Cn = P3{n} * khatri_rao(B{n}, A{n}) * pinv((B{n}'*B{n}).*(A{n}'*A{n}));
            C{n} = Stochasticize(alpha * Cn + (1-alpha) * (B{n} * T));
        end


        Ts = cell(N,1);

        % Calculate estimate of T1 as a permutation reference for other Ti's
        Ts{1} = Stochasticize(pinv(B{1}) * C{1});
        T = Ts{1}/N;

        for n = 2:N
            % Calculate Ti and permutation match to T1
            Ts{n} = Stochasticize(pinv(B{n}) * C{n});
            Pbest = PermutationFit(Ts{1}, Ts{n}, true);
            Ts{n} = Pbest * Ts{n} * Pbest';
            
            % Calculate partial T by averaging all Ti's
            T = T + Ts{n}/N;
            
            % Apply found permutation to all factor matrices
            A{n} = A{n} * Pbest';
            B{n} = B{n} * Pbest';
            C{n} = C{n} * Pbest';
        end
        
        % Update T whilst taking into consideration T_prev to reduce
        % discontinuities.
        T = rho * T_prev + (1-rho) * T;
        

        % If T is ergodic, calculate pi based on T
        pi = findStationaryDistribution(T')';

        % Set outputs
        output.tconv(iter) = norm(T-T_prev, "fro") / norm(T, "fro");
        for n = 1:N
            genPn = cpdgen({A{n}, B{n}, C{n}});
            output.recerr(n, iter) = norm(P{n} - genPn / sum(genPn, 'all'), "fro");
            output.rellrecerr(n, iter) = norm(P{n} - genPn / sum(genPn, 'all'), "fro")/norm(P{n}, "fro");
            output.emiscorr(n, iter) = sum(corrcoef(B{n}), "all");
        end


        
        % Check early stopping condition
        if(output.tconv(iter) < options.tol)
            break;
        end
    end

    output.Ts = Ts; % Output all Ti's from the last iteration
    output.tconv = output.tconv(1:iter);
    output.recerr = output.recerr(:, 1:iter);
    output.rellrecerr = output.rellrecerr(:, 1:iter);
    output.emiscorr = output.emiscorr(:, 1:iter);
end