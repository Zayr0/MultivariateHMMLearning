function A = Stochasticize(A)
A = abs(A) ./ sum(abs(A), 1);
end