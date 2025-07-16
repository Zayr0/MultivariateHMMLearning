function Z = mode_n_product(X,Y,n)
    % MODE_N_PRODUCT takes tensor X and compatible matrix Y and performs mode-n product between X and Y.
    % INPUT tensor X, matrix Y.
    % OUTPUT tensor Z.
  
    N = size(X);
    X = mode_n_matricization(X, n);
    Z = Y * X;
    Z = reshape(Z, [size(Z, 1), N(1:n-1), N(n+1:end)]);
    Z = permute(Z, [2:n, 1, n+1:length(N)]);
end