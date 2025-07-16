function P3_inf = calculateJPT(T, O, pi)
D{1} = Stochasticize(O * diag(pi) * T' / inv(diag(T * pi')));
D{2} = O;
D{3} = Stochasticize(O * T);

P3_inf = cpdgen(D);
end