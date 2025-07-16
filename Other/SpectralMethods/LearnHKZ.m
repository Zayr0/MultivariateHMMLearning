% [Hsu et al. 2012] Daniel Hsu, Sham M. Kakade and Tong Zhang. A spectral
% algorithm for learning Hidden Markov Models. Journal of Computer
% and System Sciences, vol. 78, no. 5, pages 1460–1480, September
% 2012. (Cited on pages v, 1, 2, 3, 11, 12, 13, 18, 29, 36 and 39.)



function [O_hat, T_hat, Pi_hat] = LearnHKZ(P1, P21, P31, P321, k, d)
% Data: N triples of observations, k - number of states , d - number of observations
% Result: Hidden Markov model parameterized by O_hat, T_hat and Pi_hat

[u, ~, ~] = svd(P21);
U = u(:, 1:k);

O_hat = zeros(d, k);

for r = 1:d
    eigs = eig((U' * squeeze(P321(:, r, :))) * pinv(U' * P31));
    O_hat(r, :) = real(eigs)';
end


Oinv = pinv(O_hat);

Pi_hat = Oinv*P1;


T_hat = Oinv * P21 * (Oinv') / diag(Pi_hat);

sum(O_hat, 1)

% O_hat = O_hat ./ sum(O_hat, 1);
% Pi_hat = Pi_hat ./ sum(Pi_hat);
% T_hat = T_hat ./ sum(T_hat, 1);


end