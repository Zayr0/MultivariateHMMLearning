function [T_hat, O_hat, pi_hat, Pbest] = remove_hmm_permutation_emission(O, T_hat, O_hat, pi_hat)
    Pbest = PermutationFitEmission(O, O_hat);
    T_hat = Pbest * T_hat * Pbest';
    O_hat = O_hat * Pbest';
    pi_hat = pi_hat * Pbest';
end

