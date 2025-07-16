function [M, output] = ProbabilityTensorGenerate(seq, wSize, nBins)
dims = nBins * ones(1, wSize);
M = zeros(dims);
[Y,~] = discretize(seq, nBins);

if(wSize == 3)
    for i = 1:length(Y)-2
        a = Y(i);
        b = Y(i+1);
        c = Y(i+2);
        M(a, b, c) = M(a, b, c) + 1;
    end
    output.WeightMat = M;
    M = M / sum(M(:));
else
    
    for t = 1:(length(seq) - (wSize-1))
        tensorIndex = num2cell(Y(t:t+wSize-1));
        M(tensorIndex{:}) = M(tensorIndex{:}) + 1;
    end
    
    M = M/(length(seq) - (wSize-1));

end
end