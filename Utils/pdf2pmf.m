function pmfs = pdf2pmf(O, bin_centers)
    [n, k] = size(O);
    d = length(bin_centers);
    pmfs = cell(n, 1);

    for i = 1:n
        O_hist = zeros(d, k);
        for j = 1:k
           O_hist(:, j) = Stochasticize(pdf(O{i,j}, bin_centers)')';
        end
        pmfs{i} = O_hist; 
    end
end