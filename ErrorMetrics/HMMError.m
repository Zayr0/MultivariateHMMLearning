function J = HMMError(T, O, T_hat, O_hat, plot)
    J = norm(T-T_hat, "fro") / norm(T, "fro");
    
    if(~iscell(O) && ~iscell(O_hat))
        J = J + norm(O-O_hat, "fro") / norm(O, "fro");

    elseif(iscell(O) && iscell(O_hat))
        n = length(O);
        for i = 1:n
             J = J + (1/n) * norm(O{i} - O_hat{i}, "fro") / norm(O{i}, "fro");
        end
        
    else
        warning("O and O_hat do not have the same n");
    end

    if plot
        plotHMMError(T, O, T_hat, O_hat);
    end
end