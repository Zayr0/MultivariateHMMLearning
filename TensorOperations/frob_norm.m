function z = frob_norm(X)
    % FROB_NORM takes as input tensor X of any order and returns 
    % its Frobenius norm.
    % INPUT tensor X.
    % OUTPUT scalar z.
    
    inner = frob_inner(X, X);
    z = sqrt(inner);
end