function [T_hat, O_hat, pi_hat] = FactorMatrices2HMM(U_hat)
B2 = Stochasticize(U_hat{2});
B3 = Stochasticize(U_hat{3});
T_hat = Stochasticize(pinv(B2) * B3);
O_hat = Stochasticize(B2);
pi_hat = findStationaryDistribution(T_hat');
end