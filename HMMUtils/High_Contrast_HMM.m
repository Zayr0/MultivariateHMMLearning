function [T, E] = High_Contrast_HMM(N, M, c)
    % Generates an HMM with high-contrast transition (T) and emission (E) matrices
    % N: Number of states
    % M: Number of observations
    % c: Contrast constant

    % Generate high-contrast transition matrix
    T = rand(N, N);  % Start with a random matrix
    T = T .^ c;      % Apply power to enhance contrast (makes some values very small)
    T = T ./ sum(T, 1);  % Normalize rows to sum to 1

    % Generate high-contrast emission matrix
    E = rand(M, N);  % Random values
    E = E .^ c;      % Apply power to exaggerate contrast
    E = E ./ sum(E, 1);  % Normalize rows

    % Optional: Ensure some dominant values per row
    for i = 1:N
        [~, maxIdx] = max(E(i, :)); % Find max index
        E(i, :) = E(i, :) * 0.5;    % Reduce overall values
        E(i, maxIdx) = 0.8;         % Assign a strong peak
        E(i, :) = E(i, :) / sum(E(i, :)); % Re-normalize
    end

    T = T';
    E = E';
end