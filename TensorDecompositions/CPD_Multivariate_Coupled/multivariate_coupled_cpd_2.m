function [T, B, pi, output] = multivariate_coupled_cpd_2(P, K, D, N, alpha, rho, options)
    
    % Initialize the factor matrices and transition matrix
    
    P1 = cell(N, 1); P2 = cell(N, 1); P3 = cell(N, 1);

    for n = 1:N
        P1{n} = mode_n_matricization(P{n}, 1);
        P2{n} = mode_n_matricization(P{n}, 2);
        P3{n} = mode_n_matricization(P{n}, 3);
    end
    
    if isfield(options, 'epsilon')
        epsilon = options.epsilon;
    else
        epsilon = 1e-3;
    end
    
    T = Stochasticize(eye(K) + epsilon * rand(K, K));
    pi = findStationaryDistribution(T')';

    if ~isfield(options, 'Initialfactor')   
        A = cell(N, 1); B = cell(N, 1); C = cell(N, 1);
        for n = 1:N
            A{n} = Stochasticize(ones(D,K) + epsilon * rand(D, K));
            B{n} = Stochasticize(ones(D,K) + epsilon * rand(D, K));
            C{n} = Stochasticize(ones(D,K) + epsilon * rand(D, K));
        end
    else            
        A = options.Initialfactor{1};
        B = options.Initialfactor{2};
        C = options.Initialfactor{3};
    end

    
    
    % Allocate memory for error metrics: convergence of transition matrix,
    % and tensor reconstruction error.
    output.tconv = ones(1, options.maxIter);
    output.recerr = ones(N, options.maxIter);
    output.rellrecerr = ones(N, options.maxIter);
    output.emiscorr = ones(N, options.maxIter);
    output.Ps = cell(N, 1);

    % Options for when the true variables (T, Os, Ps) are known
    output.Trec = ones(1, options.maxIter);
    output.Orec = ones(N, options.maxIter);
    output.Prec = ones(N, options.maxIter);


    for iter = 1:options.maxIter
        T_prev = T; % Save previous T to check convergence later
        
        % Calculate factor matrices of every joint probability tensor using
        % ALS with partial dependece on common T
        for n = 1:N
            A{n} = P1{n} * khatri_rao(C{n}, B{n}) * pinv((C{n}'*C{n}).*(B{n}'*B{n})) + alpha * (B{n} * diag(pi) * T' / (diag(T * pi)));
            A{n} = Stochasticize(A{n});
            
            B{n} = P2{n} * khatri_rao(C{n}, A{n}) * pinv((C{n}'*C{n}).*(A{n}'*A{n})) + alpha * B{n};
            B{n} = Stochasticize(B{n});

            C{n} = P3{n} * khatri_rao(B{n}, A{n}) * pinv((B{n}'*B{n}).*(A{n}'*A{n})) + alpha * (B{n} * T);
            C{n} = Stochasticize(C{n});
        end


        Ts = cell(N,1);

        % Calculate estimate of T1 as a permutation reference for other Ti's
        if(isfield(options, 'True_T'))
            Ts{1} = Stochasticize(pinv(B{1}) * C{1});
            Pbest = PermutationFit(options.True_T, Ts{1}, true);
            Ts{1} = Pbest * Ts{1} * Pbest';
            T = Ts{1}/N;

            % Apply found permutation to all factor matrices
            A{1} = A{1} * Pbest';
            B{1} = B{1} * Pbest';
            C{1} = C{1} * Pbest';
        else
            Ts{1} = Stochasticize(pinv(B{1}) * C{1});
            T = Ts{1}/N;
        end

        for n = 2:N
            % Calculate Ti and permutation match to T1
            Ts{n} = Stochasticize(pinv(B{n}) * C{n});

            if(isfield(options, 'True_T'))
                Pbest = PermutationFit(options.True_T, Ts{n}, true);
            else
                Pbest = PermutationFit(Ts{1}, Ts{n}, true);
            end
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

        if(isfield(options, 'True_T'))
            output.Trec(iter) = norm(options.True_T - T, "fro");
        end

        for n = 1:N
            genPn = cpdgen({A{n}, B{n}, C{n}});
            genPn = genPn ./ sum(genPn, "all");

            output.recerr(n, iter) = norm(P{n} - genPn, "fro");
            output.rellrecerr(n, iter) = norm(P{n} - genPn, "fro")/norm(P{n}, "fro");
            output.emiscorr(n, iter) = sum(corrcoef(B{n}), "all");
            
            if(isfield(options, 'True_Os'))
                output.Orec(n, iter) = norm(options.True_Os{n} - B{n}, "fro");
            end
            if(isfield(options, 'True_Ps'))
                output.Prec(n, iter) = norm(options.True_Ps{n} - genPn, "fro");
            end
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
    output.A = A;
    output.B = B;
    output.C = C;

    output.Trec = output.Trec(:, 1:iter);
    output.Orec = output.Orec(:, 1:iter);
    output.Prec = output.Prec(:, 1:iter);

    for n = 1:N
        genPn = cpdgen({A{n}, B{n}, C{n}});
        output.Ps{n} = genPn ./ sum(genPn, "all"); 
    end
end