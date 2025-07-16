function [T, O, pi, k, d] = generateContinuousHMM(k, d, n, options)

T = Stochasticize(eye(k) + options.epsilon * abs(randn(k, k)));
pi = findStationaryDistribution(T');
O = cell(n, k);

if (options.distribution == "Normal")
    sigma = 1.0;
    spacing = 2*(1:k);
    for i = 1:n
        spacing = spacing(randperm(k));
        for j = 1:k
            O{i, j} = makedist('Normal', 'mu', spacing(j), 'sigma', sigma + 0.25 * sigma * rand(1, 1));
        end
    end
end


if (options.distribution == "Gamma")
    
    for i = 1:n
        for j = 1:k
            a = 5 * rand(1) + 1;
            b = 6 * rand(1) + 4;
            O{i, j} = makedist('Gamma', 'a', a, 'b',b);
        end
    end
end




end

