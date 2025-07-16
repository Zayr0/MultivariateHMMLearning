function [G, rel_error] = TT_SVD(tensor, r)

C = tensor;
d = ndims(C);

G = cell(d);

iter = 1:d-1;
for k=iter
    C = reshape(C, [r(k-1) * n(k), numel(C)/(r(k-1) * n(k))]);

    [U, S, VT] = svd(C);
    U = U(1:r(k), 1:r(k));
    G{k} = reshape(U, [r(k-1), n(k), r(k)]);
    C = S * VT;
end
G{d} = C;
end