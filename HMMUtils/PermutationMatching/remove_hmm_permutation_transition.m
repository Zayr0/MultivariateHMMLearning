function [T_hat, O_hat, pi_hat, Pbest] = remove_hmm_permutation_transition(T, T_hat, O_hat, pi_hat)
    Pbest = PermutationFit(T, T_hat, true);
    T_hat = Pbest * T_hat * Pbest';
    pi_hat = Pbest * pi_hat;

    if(~iscell(O_hat))
        O_hat = O_hat * Pbest';
    else
        n = length(O_hat);

        for i = 1:n
            O_hat{i} = O_hat{i} * Pbest';
        end
    end
    
end

