function plotHMMError(T, O, T_hat, O_hat)

pi = findStationaryDistribution(T')';
pi_hat = findStationaryDistribution(T_hat')';

if(~iscell(O) && ~iscell(O_hat))
    plotUnivariateHMMError(T, O, pi, T_hat, O_hat, pi_hat);
end

if(iscell(O) && iscell(O_hat))
    plotMultivariateHMMError(T, O, pi, T_hat, O_hat, pi_hat);
end

if((~iscell(O) && iscell(O_hat))||(iscell(O) && ~iscell(O_hat)))
    warning("O or O_hat has the wrong variable type, check n")
end

end