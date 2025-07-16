function [Pbest, perm] = initial_permutation_fit(T, T_hat)
dT = diag(T);
dT_hat = diag(T_hat);
k = length(dT);

perm = zeros(1, k);

for i = 1:k
    [M, I] = min(dT_hat - dT(i));
    perm(i) = I;
end



if(length(perm) ~= length(unique(perm)))
    perm = NaN;
    Pbest = NaN;
else
    Pbest = eye(k);
    Pbest = permute(Pbest, perm);
end
end