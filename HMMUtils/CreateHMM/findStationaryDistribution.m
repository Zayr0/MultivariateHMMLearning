function pi = findStationaryDistribution(T)
    % Ensure T is a square matrix
    [n, m] = size(T);
    if n ~= m
        error('Transition matrix T must be square');
    end

    % Solve (T' - I) * pi' = 0 with sum(pi) = 1
    A = [T' - eye(n); ones(1, n)];
    b = [zeros(n, 1); 1];

    % Solve the system
    pi = A \ b;

    % Ensure the result is a row vector
    pi = pi';
end