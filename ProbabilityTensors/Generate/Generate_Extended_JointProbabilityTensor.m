function [M] = Generate_Extended_JointProbabilityTensor(seq, d, n)
N = 2 * n + 1;
[Y, ~] = discretize(seq, d);

M = zeros(d^n, d^n, d);

for i = (n+1):length(seq)-n
    lm1mn = Y(i-n:i-1);
    l1n = Y(i+1:i+n);
    
    l0 = Y(i);

    M(L(lm1mn, d), L(l1n, d), l0) = M(L(lm1mn, d), L(l1n, d), l0) + 1;
end

M = M / sum(M, "all");

end

function dn = L(l1n, d)
    n = length(l1n);
    dn = 0;

    for i = 1:n-1
        dn = dn + (l1n(i) - 1) * d^(n-i);
    end
    dn = dn + l1n(n);
end