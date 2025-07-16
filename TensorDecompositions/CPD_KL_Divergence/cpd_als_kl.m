function [B1, B2, B3, T, O] = cpd_als_kl(P, k, lambda_T, lambda_O, max_iter, eta)
% ALS_KL_HMM: Tensor Decomposition for HMMs with KL Regularization
%   P: dxdxd joint probability tensor
%   k: Number of hidden states
%   lambda_T: KL regularization weight for T
%   lambda_O: KL regularization weight for O
%   max_iter: Number of ALS iterations
%   eta: Learning rate for KL updates
%   Returns: Estimated Transition (T) and Emission (O) matrices

[d, ~, ~] = size(P);

% Initialize factor matrices using SVD
[B1, B2, B3] = initialize_factors(P, k);

% Initialize T and O from factor matrices
T = Stochasticize(rand(k, k));
O = Stochasticize(rand(d, k));

% Define priors for KL regularization
T_prior = Stochasticize(eye(k) + 0.1 * rand(k, k));
O_prior = Stochasticize(ones(d, k) / d);

P1 = unfold(P, 1);
P2 = unfold(P, 2);
P3 = unfold(P, 3);

for iter = 1:max_iter

    B2B3_kr = khatri_rao(B3, B2);
    B1B3_kr = khatri_rao(B3, B1);
    B1B2_kr = khatri_rao(B2, B1);


    % ALS Updates
    B1 = P1 * B2B3_kr / (B2B3_kr' * B2B3_kr);
    B1 = Stochasticize(B1);

    B2 = P2 * B1B3_kr / (B1B3_kr' * B1B3_kr);
    B2 = Stochasticize(B2);

    B3 = P3 * B1B2_kr / (B1B2_kr' * B1B2_kr);
    B3 = Stochasticize(B3);


    T = pinv(B2'*B2 + lambda_T * eye(k)) * (B2' * B3 + lambda_T * T_prior);
    O = (B2 + lambda_O * O_prior)/(1+lambda_O);
    
    % % Extract T and O from factor matrices
    % T = pinv(B2) * B3;
    % O = B2;
    % 
    % % KL Regularization Updates
    % T = T - eta * lambda_T * (log(T ./ T_prior) + 1);
    % O = O - eta * lambda_O * (log(O ./ O_prior) + 1);
    % 
    % % Ensure column stochasticity
    T = Stochasticize(T);
    O = Stochasticize(O);
end

end

function B_new = update_factor(P_i, B1, B2)
    B_new = P_i * khatri_rao(B2, B1) / ((B2' * B2) .* (B1' * B1));
end

function X = normalize_columns(X)
    X = X ./ sum(X, 1);
end

function [B1, B2, B3] = initialize_factors(P, k)
    [U, S, V] = svd(unfold(P, 1), 'econ');
    B1 = U(:, 1:k) * sqrt(S(1:k, 1:k));
    B2 = (V(1:k, :)') * sqrt(S(1:k, 1:k));
    B3 = rand(size(B2));
end

% Helper function to unfold tensor along mode-n
function Xn = unfold(X, mode)
    sz = size(X);
    Xn = reshape(permute(X, [mode, setdiff(1:ndims(X), mode)]), sz(mode), []);
end
