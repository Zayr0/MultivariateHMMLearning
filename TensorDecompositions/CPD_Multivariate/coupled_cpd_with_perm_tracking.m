function [A, B, C] = coupled_cpd_with_perm_tracking(P_tensors, k, max_iter)
    % P_tensors: Cell array of joint probability tensors {P1, P2, ..., Pn}
    % k: Number of latent states
    % max_iter: Maximum ALS iterations
    % A, B, C: Factor matrices of size d x (n*k)
    % perm: Permutation tracking matrix
    
    n = length(P_tensors);  % Number of time series
    d = size(P_tensors{1}, 1);  % Emission dimension
    
    % Initialize factor matrices with random values
    A = rand(d, n*k);
    B = rand(d, n*k);
    C = rand(d, n*k);

    
    for iter = 1:max_iter
        % Update A
        for i = 1:n
            Pi = P_tensors{i};
            unfold_Pi = reshape(Pi, d, []);  % Mode-1 unfolding
            B_i = B(:, (i-1)*k+1:i*k);
            C_i = C(:, (i-1)*k+1:i*k);
            
            kr_BC = khatri_rao(C_i, B_i);  % (d^2 x k)
            denom = (C_i' * C_i) .* (B_i' * B_i); % k x k
            
            % Update A and track permutation
            A(:, (i-1)*k+1:i*k) = ((unfold_Pi * kr_BC) / denom);
        end
        
        % Update B
        for i = 1:n
            Pi = P_tensors{i};
            unfold_Pi = reshape(permute(Pi, [2, 1, 3]), d, []);  % Mode-2 unfolding
            A_i = A(:, (i-1)*k+1:i*k);
            C_i = C(:, (i-1)*k+1:i*k);
            
            kr_AC = khatri_rao(C_i, A_i);
            denom = (C_i' * C_i) .* (A_i' * A_i);
            
            % Update B and track permutation
            B(:, (i-1)*k+1:i*k) = ((unfold_Pi * kr_AC) / denom);
        end
        
        % Update C
        for i = 1:n
            Pi = P_tensors{i};
            unfold_Pi = reshape(permute(Pi, [3, 1, 2]), d, []);  % Mode-3 unfolding
            A_i = A(:, (i-1)*k+1:i*k);
            B_i = B(:, (i-1)*k+1:i*k);
            
            kr_AB = khatri_rao(B_i, A_i);
            denom = (A_i' * A_i) .* (B_i' * B_i);
            
            % Update C and track permutation
            C(:, (i-1)*k+1:i*k) = ((unfold_Pi * kr_AB) / denom);
        end
    end
end
