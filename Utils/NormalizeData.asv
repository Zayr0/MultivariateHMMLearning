function [seq, mins, maxs] = NormalizeData(seq)
% The aim of this function is to scale all data between zero and one.
% Each data sequence has its own range, this range needs to be saved to
% backtrack the orignal values

% seq is an n x m data matrix where n is the number of sequences and m is
% the length of each sequence.
    [n, m] = size(seq);

    mins = min(seq, [], 2);
    maxs = max(seq, [], 2);

    seq = 2 * (seq - mins * ones(1, m)) ./ ((maxs - mins) * ones(1, m)) - 1;
    
end