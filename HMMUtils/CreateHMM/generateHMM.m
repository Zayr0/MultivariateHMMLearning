function [T, O, pi, k, d] = generateHMM(k, d, n, options)

if(options.method == "random")
    T = eye(k) + options.epsilon * rand(k, k);
    T = Stochasticize(T);
    O = initialize_emission_matrix_offset(d, k, n);

elseif(options.method == "random_no_offset")
    T = eye(k) + options.epsilon * rand(k, k);
    T = Stochasticize(T);
    O = initialize_emission_matrix(d, k, n);
    
elseif((options.method == "standard") && n == 1)
    k = 3;
    d = 8;
    
    T = [8/10, 1/15, 1/6; 
        1/10, 13/15, 1/6; 
        1/10, 1/15, 2/3];
    
    O = [6/15, 1/20, 1/50;
        1/15, 11/20, 1/50;
        1/15, 1/20, 1/50;
        1/15, 1/20, 21/50;
        1/15, 1/20, 21/50;
        1/15, 1/20, 1/50;
        1/15, 1/20, 1/50;
        3/15, 3/20, 3/50;];
elseif((options.method == "standard") && n ~= 1)
    warning("Standard HMM cannot generate multiple observational matrices.")
end
pi = findStationaryDistribution(T');
end

