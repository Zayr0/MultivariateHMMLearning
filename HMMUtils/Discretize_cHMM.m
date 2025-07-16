function [T, O_hists] = Discretize_cHMM(T, Os, mins, maxs, d)
    n = length(mins);
    k = size(T, 1);
    
    t_est = zeros(n, d);
    range = maxs - mins;
    O_hists = cell(n, 1);
    
    for i = 1:n
        t_est(i, :) = linspace(mins(i), maxs(i), d);
        O_hist = zeros(d, k);
        
        for j = 1:k
            y = pdf(Os{i,j}, t_est(i, :)) * range(i) / d;
            O_hist(:, j) = y;
        end
        O_hists{i} = O_hist;
    end
end

