function [A, B, C, D, conv] = JPT_NTF(T, M, k, n, options)


A = rand(n, k);
B = rand(n, k);
C = rand(n, k);
D = rand(k, k);

T1 = mode_n_matricization(T, 1);
T2 = mode_n_matricization(T, 2);
T3 = mode_n_matricization(T, 3);

conv = zeros(1, options.maxIter);

for i = 1:options.maxIter
    A = A .* ((T1 * khatri_rao(C, B)) ./ ((A*((C'*C).*(B'*B)))));
    B = B .* ((T2 * khatri_rao(C, A) + options.rho * (M*B*D' + M'*B*D)) ./ (B*((A'*A).*(C'*C)) + options.rho * (B*D'*(B')*B*D + B*D*(B')*B*D')));
    C = C .* ((T3 * khatri_rao(B, A)) ./ (C*((B'*B).*(A'*A))));
    D = D .* ((B'*M*B) ./ (B'*B*D*(B')*B));

    J = 0.5 * (norm(T - reshape(A * khatri_rao(C, B)', n, n, n), "fro")^2 + norm(M - B * D * B', "fro")^2);
    conv(i) = J;
    
    if(J < options.epsilon)
        display("Early NTF termination, minimum error reached at iteration: " + i);
        break;
    end

    if(i > 1)
        if(abs(conv(i) - conv(i-1)) < options.delta)
            display("Early NTF termination, minimal change between iterations ");
            break;
        end
    end
end

i = 1;

A = A * diag(sum(B, i)) * diag(sum(C, i));
A = A ./ sum(A, "all");

D = diag(sum(B, i)) * D * diag(sum(B, i));
D = D ./ sum(D, "all");

B = B * inv(diag(sum(B, i)));
C = C * inv(diag(sum(C, i)));

% A = A ./ sum(A, 1);
% B = B ./ sum(B, 1);
% C = C ./ sum(C, 1);
% D = D ./ sum(D, 1);

end