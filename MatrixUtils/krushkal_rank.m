function Kr = krushkal_rank(M)
    R = rank(M);
    Kr = R;

    for i = 1:size(M, 2)
        v = M(:, i); % Column vector to check
        S = M(:, 1:end ~= i); % Matrix M with column vector excluded

        A = sym("a",[size(S, 2) 1]); % Symbolic variables
        eqns = (sum(S*A)-v == 0); % Create symbolic equation to solving
        b = struct2cell(solve(eqns,A, Real=true)); % Solve equation
        n = length(find([b{:}]~=0)); % Find the number of vectors that v is linearly dependent on
        
        if n < Kr && n ~= 0
            Kr = n; % Kruskal rank is equal to the minimum n, for all linearly dependent vectors v_i
        end
    end
end