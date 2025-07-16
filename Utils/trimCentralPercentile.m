function trimmed_seq = trimCentralPercentile(seq, keep_percent)
    % Removes the tails, keeping only the central 'keep_percent' of the data.
    % For 99%, 0.5% is removed from each end.
    if nargin < 2
        keep_percent = 99;
    end

    lower_cutoff = (100 - keep_percent) / 2;
    upper_cutoff = 100 - lower_cutoff;

    lower_bound = prctile(seq, lower_cutoff);
    upper_bound = prctile(seq, upper_cutoff);

    trimmed_seq = seq(seq >= lower_bound & seq <= upper_bound);
end