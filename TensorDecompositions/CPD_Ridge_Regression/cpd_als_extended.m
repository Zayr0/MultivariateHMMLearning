function [B1, B2, B3, output] = cpd_als_extended(P, k, options)
    % CP Decomposition with Structured Ridge Regression ALS for HMMs
    % Inputs:
    %   P       - Joint probability tensor (d x d x d)
    %   k       - Number of states (rank of decomposition)
    %   options.lambda1 - Regularization weight for B1 ~ B2 N
    %   options.lambda2 - Regularization weight for B3 ~ B2 M
    %   options.maxIter - Maximum number of ALS iterations
    %   options.tol     - Convergence options.tolerance
    % Outputs:
    %   B1, B2, B3 - Factor matrices of the CP decomposition

    % Get tensor size
    d = size(P, 1);

    epsilon = 1e-3;
    T_init = Stochasticize(eye(k) + epsilon * rand(k,k));
    O_init = Stochasticize(rand(d, k));
    pi_init = findStationaryDistribution(T_init')';

    B1 = O_init * diag(pi_init) * T_init' * inv(diag(T_init*pi_init));
    B2 = O_init;
    B3 = O_init * T_init;

    % % Initialize factor matrices randomly
    % B1 = rand(d, k);
    % B2 = rand(d, k);
    % B3 = rand(d, k);

    P1 = mode_n_matricization(P, 1);
    P2 = mode_n_matricization(P, 2);
    P3 = mode_n_matricization(P, 3);

    D = diff(eye(d));
    DtD = D' * D;

    output.relerr_tensor = ones(1, options.maxIter);
    output.convergence_factor = ones(3, options.maxIter);
    output.cost = ones(1, options.maxIter);

    % Main ALS loop
    for iter = 1:options.maxIter
        % Store previous values to check for convergence
        B1_old = B1;
        B2_old = B2;
        B3_old = B3;
        
        % Compute pseudoinverses for structured constraints
        B2_pinv = pinv(B2);  % Pseudoinverse of B2
        
        % Compute Khatri-Rao products for ALS updates
        % B2B3_kr = khatri_rao(B3, B2);
        % B1B3_kr = khatri_rao(B3, B1);
        % B1B2_kr = khatri_rao(B2, B1);
    
        if(options.mode == "Regular")
            B1 = P1 * khatri_rao(B3, B2) * pinv((B3'*B3).*(B2'*B2));
            B1 = Stochasticize(B1);
        
            B2 = P2 * khatri_rao(B3, B1) * pinv((B3'*B3).*(B1'*B1));
            B2 = Stochasticize(B2);
            
            B3 = P3 * khatri_rao(B2, B1) * pinv((B2'*B2).*(B1'*B1));
            B3 = Stochasticize(B3);

        elseif(options.mode == "Structured")
            % Update N and M
            % N = Stochasticize((B2_pinv * B1)')';
            % M = Stochasticize(B2_pinv * B3);

            N = pinv(B2' * B2) * B2' * B1;
            M = pinv(B2' * B2) * B2' * B3;

            N = Stochasticize(N')';
            M = Stochasticize(M);

            % Update B1
            B1 = (P1 * khatri_rao(B3, B2) - options.lambda1 * B2 * N) * pinv((B3'*B3).*(B2'*B2) + options.lambda1 * eye(k));
            B1 = Stochasticize(B1);
    
            % Update B2
            B2 = (P2 * khatri_rao(B3, B1) + options.lambda1 * B1 * N' + options.lambda2 * B3 * M') * pinv((B3'*B3).*(B1'*B1) + (options.lambda1 + options.lambda2) * eye(k));
            B2 = Stochasticize(B2);
    
            % Update B3
            B3 = (P3 * khatri_rao(B2, B1) - options.lambda2 * B2 * M) * pinv((B2'*B2).*(B1'*B1) + options.lambda2 * eye(k));
            B3 = Stochasticize(B3);

        elseif(options.mode == "Smooth")
            % Update B1
            B1 = (P1 * khatri_rao(B3, B2) + options.lambda1 * DtD * B1) *  pinv((B3'*B3).*(B2'*B2) + options.lambda1 * eye(k));
            B1 = Stochasticize(B1);
    
            % Update B2
            B2 = (P2 * khatri_rao(B3, B1) + options.lambda1 * DtD * B2) * pinv((B3'*B3).*(B1'*B1) + options.lambda1 * eye(k));
            B2 = Stochasticize(B2);
    
            % Update B3
            B3 = (P3 * khatri_rao(B2, B1) + options.lambda1 * DtD * B3) * pinv((B2'*B2).*(B1'*B1) + options.lambda1 * eye(k));
            B3 = Stochasticize(B3);

        elseif(options.mode == "Ridge")
            % Update B1
            B1 = (P1 * khatri_rao(B3, B2)) * pinv((B3'*B3).*(B2'*B2) + options.lambda1 * eye(k));
            B1 = Stochasticize(B1);
    
            % Update B2
            B2 = (P2 * khatri_rao(B3, B1)) * pinv((B3'*B3).*(B1'*B1) + options.lambda1 * eye(k));
            B2 = Stochasticize(B2);
    
            % Update B3
            B3 = (P3 * khatri_rao(B2, B1)) * pinv((B2'*B2).*(B1'*B1) + options.lambda1 * eye(k));
            B3 = Stochasticize(B3);

        elseif(options.mode == "Proximal")
            % Update B1
            B1 = (P1 * khatri_rao(B3, B2) + options.lambda1 * B1_old) * pinv((B3'*B3).*(B2'*B2) + options.lambda1 * eye(k));
            B1 = Stochasticize(B1);

            % Update B2
            B2 = (P2 * khatri_rao(B3, B1) + options.lambda1 * B2_old) * pinv((B3'*B3).*(B1'*B1) + options.lambda1  * eye(k));
            B2 = Stochasticize(B2);
    
            % Update B3
            B3 = (P3 * khatri_rao(B2, B1) + options.lambda1 * B3_old) * pinv((B2'*B2).*(B1'*B1) + options.lambda1 * eye(k));
            B3 = Stochasticize(B3);
        elseif(options.mode == "Prior")
            % Update B1
            B1 = (P1 * khatri_rao(B3, B2) - options.lambda1 * B2 * eye(k)) * pinv((B3'*B3).*(B2'*B2) + options.lambda1 * eye(k));
            B1 = Stochasticize(B1);

            % Update B2
            B2 = (P2 * khatri_rao(B3, B1) + options.lambda1 * B1 * eye(k) + options.lambda2 * B3 * eye(k)) * pinv((B3'*B3).*(B1'*B1)+  2 * options.lambda1  * eye(k));
            B2 = Stochasticize(B2);
    
            % Update B3
            B3 = (P3 * khatri_rao(B2, B1) - options.lambda1 * B2 * eye(k)) * pinv((B2'*B2).*(B1'*B1) + options.lambda1 * eye(k));
            B3 = Stochasticize(B3);
        else
            warning("Wrong mode selected.");
        end

        output.relerr_tensor(iter) = norm(P - cpdgen({B1, B2, B3}), "fro")^2;
        output.convergence_factor(:, iter) = [norm(B1 - B1_old, 'fro'); norm(B2 - B2_old, 'fro'); norm(B3 - B3_old, 'fro')];
        % output.cost(iter) = norm(P - cpdgen({B1, B2, B3}), "fro")^2 + options.lambda1 * norm(B1 - B2 * N, "fro")^2 + options.lambda2 * norm(B3 - B2 * M, "fro")^2;

        % Check convergence (Frobenius norm difference)
        if norm(B1 - B1_old, 'fro') < options.tol && ...
           norm(B2 - B2_old, 'fro') < options.tol && ...
           norm(B3 - B3_old, 'fro') < options.tol
            fprintf('Converged at iteration %d\n', iter);
            break;
        end
    end

    B1 = Stochasticize(B1);
    B2 = Stochasticize(B2);
    B3 = Stochasticize(B3);

    output.relerr_tensor = output.relerr_tensor(1:iter);
    output.convergence_factor = output.convergence_factor(:, 1:iter);
    output.iters = 1:iter;
end

% Helper function to unfold tensor along mode-n
function Xn = unfold(X, mode)
    sz = size(X);
    Xn = reshape(permute(X, [mode, setdiff(1:ndims(X), mode)]), sz(mode), []);
end

% Khatri-Rao Product (column-wise Kronecker product)
function KR = kr(A, B)
    [rA, cA] = size(A);
    [rB, cB] = size(B);
    if cA ~= cB
        error('Matrix dimensions must match for Khatri-Rao product');
    end
    KR = zeros(rA * rB, cA);
    for i = 1:cA
        KR(:, i) = kron(A(:, i), B(:, i));
    end
end
